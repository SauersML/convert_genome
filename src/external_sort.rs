use std::{
    cmp::Ordering,
    collections::BinaryHeap,
    fs::File,
    io::{self, BufReader, BufWriter, Write},
};

use noodles::vcf::variant::io::Write as VariantRecordWrite;
use noodles::{bcf, vcf};
use tempfile::NamedTempFile;

use crate::dtc;
use crate::dtc::Record as DtcRecord;
use crate::input::natural_contig_order;
use crate::reference::ReferenceGenome;

const RECORDBUF_CHUNK_SIZE: usize = 50_000;
const DTC_CHUNK_SIZE: usize = 500_000;

#[derive(Clone, Copy, Debug)]
pub enum SortFormat {
    Vcf,
    Bcf,
}

#[derive(Clone, Debug)]
pub enum RecordOrder {
    Reference(std::collections::HashMap<String, usize>),
    Natural,
}

impl RecordOrder {
    pub fn from_reference(reference: &ReferenceGenome) -> Self {
        Self::Reference(reference.contig_index_map())
    }

    fn key(&self, record: &vcf::variant::record_buf::RecordBuf) -> RecordSortKey {
        let pos = record.variant_start().map(usize::from).unwrap_or(0);
        match self {
            Self::Reference(map) => {
                let idx = map
                    .get(record.reference_sequence_name())
                    .copied()
                    .unwrap_or(usize::MAX);
                if idx == usize::MAX {
                    let (order, name) = natural_contig_order(record.reference_sequence_name());
                    RecordSortKey::Natural { order, name, pos }
                } else {
                    RecordSortKey::Reference { idx, pos }
                }
            }
            Self::Natural => {
                let (order, name) = natural_contig_order(record.reference_sequence_name());
                RecordSortKey::Natural { order, name, pos }
            }
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
enum RecordSortKey {
    Reference {
        idx: usize,
        pos: usize,
    },
    Natural {
        order: u32,
        name: String,
        pos: usize,
    },
}

impl Ord for RecordSortKey {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (Self::Reference { idx: a, pos: ap }, Self::Reference { idx: b, pos: bp }) => {
                (a, ap).cmp(&(b, bp))
            }
            (
                Self::Natural {
                    order: a,
                    name: an,
                    pos: ap,
                },
                Self::Natural {
                    order: b,
                    name: bn,
                    pos: bp,
                },
            ) => (a, an, ap).cmp(&(b, bn, bp)),
            (Self::Reference { .. }, Self::Natural { .. }) => Ordering::Less,
            (Self::Natural { .. }, Self::Reference { .. }) => Ordering::Greater,
        }
    }
}

impl PartialOrd for RecordSortKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

struct HeapItem {
    key: RecordSortKey,
    spill_index: usize,
    record: vcf::variant::record_buf::RecordBuf,
}

impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key.cmp(&other.key)
    }
}

impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for HeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key
    }
}

impl Eq for HeapItem {}

enum SpillReader {
    Vcf {
        reader: vcf::io::Reader<BufReader<File>>,
        header: vcf::Header,
    },
    Bcf {
        reader: bcf::io::Reader<Box<dyn std::io::Read>>,
        header: vcf::Header,
    },
}

impl SpillReader {
    fn open(spill: &NamedTempFile, format: SortFormat) -> io::Result<Self> {
        let file = spill.reopen()?;
        let reader = BufReader::new(file);
        match format {
            SortFormat::Vcf => {
                let mut reader = vcf::io::Reader::new(reader);
                let header = reader.read_header()?;
                Ok(Self::Vcf { reader, header })
            }
            SortFormat::Bcf => {
                let mut reader = bcf::io::reader::Builder::default().build_from_reader(reader)?;
                let header = reader.read_header()?;
                Ok(Self::Bcf { reader, header })
            }
        }
    }

    fn read_record(&mut self) -> io::Result<Option<vcf::variant::record_buf::RecordBuf>> {
        let mut record = vcf::variant::record_buf::RecordBuf::default();
        match self {
            Self::Vcf { reader, header } => match reader.read_record_buf(header, &mut record) {
                Ok(0) => Ok(None),
                Ok(_) => Ok(Some(record)),
                Err(e) => Err(e),
            },
            Self::Bcf { reader, header } => match reader.read_record_buf(header, &mut record) {
                Ok(0) => Ok(None),
                Ok(_) => Ok(Some(record)),
                Err(e) => Err(e),
            },
        }
    }
}

pub enum SortedRecordIter {
    InMemory(std::vec::IntoIter<vcf::variant::record_buf::RecordBuf>),
    External(RecordMergeIter),
}

impl Iterator for SortedRecordIter {
    type Item = io::Result<vcf::variant::record_buf::RecordBuf>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::InMemory(iter) => iter.next().map(Ok),
            Self::External(iter) => iter.next(),
        }
    }
}

pub struct RecordExternalSorter {
    header: vcf::Header,
    format: SortFormat,
    order: RecordOrder,
    buffer: Vec<vcf::variant::record_buf::RecordBuf>,
    spills: Vec<NamedTempFile>,
}

impl RecordExternalSorter {
    pub fn new(header: vcf::Header, format: SortFormat, order: RecordOrder) -> Self {
        Self {
            header,
            format,
            order,
            buffer: Vec::new(),
            spills: Vec::new(),
        }
    }

    pub fn push(&mut self, record: vcf::variant::record_buf::RecordBuf) -> io::Result<()> {
        self.buffer.push(record);
        if self.buffer.len() >= RECORDBUF_CHUNK_SIZE {
            self.spill()?;
        }
        Ok(())
    }

    pub fn finish(mut self) -> io::Result<SortedRecordIter> {
        if !self.spills.is_empty() {
            if !self.buffer.is_empty() {
                self.spill()?;
            }
            let iter = RecordMergeIter::new(self.spills, self.format, self.order)?;
            Ok(SortedRecordIter::External(iter))
        } else {
            self.buffer
                .sort_by(|a, b| self.order.key(a).cmp(&self.order.key(b)));
            Ok(SortedRecordIter::InMemory(self.buffer.into_iter()))
        }
    }

    fn spill(&mut self) -> io::Result<()> {
        self.buffer
            .sort_by(|a, b| self.order.key(a).cmp(&self.order.key(b)));

        let mut spill = NamedTempFile::new()?;
        match self.format {
            SortFormat::Vcf => {
                let mut writer = vcf::io::Writer::new(BufWriter::new(spill.as_file_mut()));
                writer.write_header(&self.header)?;
                for record in &self.buffer {
                    writer.write_variant_record(&self.header, record)?;
                }
            }
            SortFormat::Bcf => {
                let mut writer =
                    bcf::io::writer::Builder::default().build_from_writer(spill.as_file_mut());
                writer.write_header(&self.header)?;
                for record in &self.buffer {
                    writer.write_variant_record(&self.header, record)?;
                }
            }
        }

        self.spills.push(spill);
        self.buffer.clear();
        Ok(())
    }
}

pub struct RecordMergeIter {
    readers: Vec<SpillReader>,
    heap: BinaryHeap<std::cmp::Reverse<HeapItem>>,
    order: RecordOrder,
}

impl RecordMergeIter {
    fn new(spills: Vec<NamedTempFile>, format: SortFormat, order: RecordOrder) -> io::Result<Self> {
        let mut readers = Vec::with_capacity(spills.len());
        for spill in &spills {
            readers.push(SpillReader::open(spill, format)?);
        }

        let mut heap = BinaryHeap::new();
        for (index, reader) in readers.iter_mut().enumerate() {
            if let Some(record) = reader.read_record()? {
                let key = order.key(&record);
                heap.push(std::cmp::Reverse(HeapItem {
                    key,
                    spill_index: index,
                    record,
                }));
            }
        }

        Ok(Self {
            readers,
            heap,
            order,
        })
    }
}

impl Iterator for RecordMergeIter {
    type Item = io::Result<vcf::variant::record_buf::RecordBuf>;

    fn next(&mut self) -> Option<Self::Item> {
        let std::cmp::Reverse(item) = self.heap.pop()?;
        let spill_index = item.spill_index;
        let record = item.record;

        if let Some(reader) = self.readers.get_mut(spill_index) {
            match reader.read_record() {
                Ok(Some(next_record)) => {
                    let key = self.order.key(&next_record);
                    self.heap.push(std::cmp::Reverse(HeapItem {
                        key,
                        spill_index,
                        record: next_record,
                    }));
                }
                Ok(None) => {}
                Err(e) => return Some(Err(e)),
            }
        }

        Some(Ok(record))
    }
}

#[derive(Clone, Debug)]
pub enum DtcOrder {
    Reference(std::collections::HashMap<String, usize>),
    Natural,
}

impl DtcOrder {
    pub fn from_reference(reference: &ReferenceGenome) -> Self {
        Self::Reference(reference.contig_index_map())
    }

    fn key(&self, record: &DtcRecord) -> DtcSortKey {
        match self {
            Self::Reference(map) => {
                let idx = map.get(&record.chromosome).copied().unwrap_or(usize::MAX);
                DtcSortKey::Reference {
                    idx,
                    pos: record.position,
                }
            }
            Self::Natural => {
                let (order, name) = natural_contig_order(&record.chromosome);
                DtcSortKey::Natural {
                    order,
                    name,
                    pos: record.position,
                }
            }
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
enum DtcSortKey {
    Reference { idx: usize, pos: u64 },
    Natural { order: u32, name: String, pos: u64 },
}

impl Ord for DtcSortKey {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (Self::Reference { idx: a, pos: ap }, Self::Reference { idx: b, pos: bp }) => {
                (a, ap).cmp(&(b, bp))
            }
            (
                Self::Natural {
                    order: a,
                    name: an,
                    pos: ap,
                },
                Self::Natural {
                    order: b,
                    name: bn,
                    pos: bp,
                },
            ) => (a, an, ap).cmp(&(b, bn, bp)),
            (Self::Reference { .. }, Self::Natural { .. }) => Ordering::Less,
            (Self::Natural { .. }, Self::Reference { .. }) => Ordering::Greater,
        }
    }
}

impl PartialOrd for DtcSortKey {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

struct DtcHeapItem {
    key: DtcSortKey,
    spill_index: usize,
    record: DtcRecord,
}

impl Ord for DtcHeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key.cmp(&other.key)
    }
}

impl PartialOrd for DtcHeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for DtcHeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key
    }
}

impl Eq for DtcHeapItem {}

pub enum DtcSortedIter {
    InMemory(std::vec::IntoIter<DtcRecord>),
    External(DtcMergeIter),
}

impl Iterator for DtcSortedIter {
    type Item = DtcRecord;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::InMemory(iter) => iter.next(),
            Self::External(iter) => iter.next(),
        }
    }
}

pub struct DtcExternalSorter {
    order: DtcOrder,
    buffer: Vec<DtcRecord>,
    spills: Vec<NamedTempFile>,
}

impl DtcExternalSorter {
    pub fn new(order: DtcOrder) -> Self {
        Self {
            order,
            buffer: Vec::new(),
            spills: Vec::new(),
        }
    }

    pub fn push(&mut self, record: DtcRecord) -> io::Result<()> {
        self.buffer.push(record);
        if self.buffer.len() >= DTC_CHUNK_SIZE {
            self.spill()?;
        }
        Ok(())
    }

    pub fn finish(mut self) -> io::Result<DtcSortedIter> {
        if !self.spills.is_empty() {
            if !self.buffer.is_empty() {
                self.spill()?;
            }
            let iter = DtcMergeIter::new(self.spills, self.order)?;
            Ok(DtcSortedIter::External(iter))
        } else {
            self.buffer
                .sort_by(|a, b| self.order.key(a).cmp(&self.order.key(b)));
            Ok(DtcSortedIter::InMemory(self.buffer.into_iter()))
        }
    }

    fn spill(&mut self) -> io::Result<()> {
        self.buffer
            .sort_by(|a, b| self.order.key(a).cmp(&self.order.key(b)));

        let mut spill = NamedTempFile::new()?;
        {
            let mut writer = BufWriter::new(spill.as_file_mut());
            for record in &self.buffer {
                writeln!(writer, "{}", record)?;
            }
        }

        self.spills.push(spill);
        self.buffer.clear();
        Ok(())
    }
}

pub struct DtcMergeIter {
    readers: Vec<dtc::Reader<BufReader<File>>>,
    heap: BinaryHeap<std::cmp::Reverse<DtcHeapItem>>,
    order: DtcOrder,
}

impl DtcMergeIter {
    fn new(spills: Vec<NamedTempFile>, order: DtcOrder) -> io::Result<Self> {
        let mut readers = Vec::with_capacity(spills.len());
        for spill in &spills {
            let file = spill.reopen()?;
            let reader = BufReader::new(file);
            readers.push(dtc::Reader::new(reader));
        }

        let mut heap = BinaryHeap::new();
        for (index, reader) in readers.iter_mut().enumerate() {
            if let Some(record) = read_next_dtc(reader) {
                let key = order.key(&record);
                heap.push(std::cmp::Reverse(DtcHeapItem {
                    key,
                    spill_index: index,
                    record,
                }));
            }
        }

        Ok(Self {
            readers,
            heap,
            order,
        })
    }
}

impl Iterator for DtcMergeIter {
    type Item = DtcRecord;

    fn next(&mut self) -> Option<Self::Item> {
        let std::cmp::Reverse(item) = self.heap.pop()?;
        let spill_index = item.spill_index;
        let record = item.record;

        if let Some(reader) = self.readers.get_mut(spill_index) {
            if let Some(next_record) = read_next_dtc(reader) {
                let key = self.order.key(&next_record);
                self.heap.push(std::cmp::Reverse(DtcHeapItem {
                    key,
                    spill_index,
                    record: next_record,
                }));
            }
        }

        Some(record)
    }
}

fn read_next_dtc(reader: &mut dtc::Reader<BufReader<File>>) -> Option<DtcRecord> {
    for result in reader {
        match result {
            Ok(record) => return Some(record),
            Err(e) => {
                tracing::warn!("failed to read spill DTC record: {}", e);
            }
        }
    }
    None
}
