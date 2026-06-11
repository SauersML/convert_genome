use std::{
    fmt,
    io::{self, BufRead},
    num::ParseIntError,
};

use thiserror::Error;

/// Optional per-sample array metrics carried alongside a genotype.
///
/// Populated from an Illumina GenomeStudio Final Report (`Log R Ratio`,
/// `B Allele Freq`, `GC Score`, `GT Score`); `None` for plain DTC text inputs
/// that carry no such signal. These flow through to per-sample VCF FORMAT
/// fields so the intensity/quality information the array produced is not lost.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ArrayMetrics {
    /// B Allele Frequency (0..=1), normalized allelic intensity ratio.
    pub baf: Option<f32>,
    /// Log R Ratio, the copy-number intensity signal.
    pub lrr: Option<f32>,
    /// Illumina GenCall confidence score (the `GC Score` column, 0..=1).
    pub gencall: Option<f32>,
    /// Illumina GenTrain cluster-quality score (the `GT Score` column).
    pub gentrain: Option<f32>,
}

impl ArrayMetrics {
    /// True when every metric is absent (nothing worth carrying).
    pub fn is_empty(&self) -> bool {
        self.baf.is_none()
            && self.lrr.is_none()
            && self.gencall.is_none()
            && self.gentrain.is_none()
    }
}

/// A single genotype entry from a direct-to-consumer text export or a
/// GenomeStudio Final Report.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    pub id: Option<String>,
    pub chromosome: String,
    pub position: u64,
    pub genotype: String,
    /// Optional array metrics (GenomeStudio); `None` for DTC text inputs.
    pub metrics: Option<ArrayMetrics>,
}

impl Record {
    pub fn is_missing(&self) -> bool {
        self.genotype.trim().is_empty() || self.genotype == "--"
    }

    pub fn parse_alleles(&self) -> Result<Vec<Allele>, crate::conversion::RecordConversionError> {
        // We return Vec<Allele>. Error type?
        // Logic currently doesn't fail, just returns Missing for garbage.
        // So Ok(...) always?
        Ok(parse_genotype(&self.genotype))
    }
}

/// Represents an allele from DTC genotype data.
/// DTC data can contain SNPs (A,C,G,T), deletions (D), insertions (I), or missing (-).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Allele {
    Base(String),
    Deletion,
    Insertion,
    Missing,
}

/// Parse a genotype string into Allele states.
pub fn parse_genotype(raw: &str) -> Vec<Allele> {
    let trimmed = raw.trim();

    if trimmed.contains('/') {
        trimmed
            .split('/')
            .map(|s| match s.trim().to_ascii_uppercase().as_str() {
                "D" => Allele::Deletion,
                "I" => Allele::Insertion,
                "-" | "0" | "?" | "" => Allele::Missing,
                // Uppercase the base: downstream genotype formatting compares
                // case-sensitively against the uppercase REF/ALT, so a lowercase
                // call (e.g. "a/g") would otherwise be dropped as invalid.
                val => Allele::Base(val.to_string()),
            })
            .collect()
    } else {
        // No separator.
        // Heuristic: If string is > 1 char and contains ONLY ACGTN, treat as single Indel allele.
        // If it contains I, D, -, or is 2 chars like "AG" (diploid SNV), check against standard known calls.

        let upper = trimmed.to_ascii_uppercase();
        let is_all_bases = upper
            .chars()
            .all(|c| matches!(c, 'A' | 'C' | 'G' | 'T' | 'N'));

        // Treat as single string (Insertion/Indel) only if length > 2
        // Length 2 (e.g. "AA", "AG") is typically a diploid call.
        // Length 3+ (e.g. "ACT") cannot be diploid (triploid?), so it must be an insertion sequence.
        if upper.len() > 2 && is_all_bases {
            vec![Allele::Base(upper)]
        } else {
            // Standard splitting (e.g. "AG" -> A, G; "II" -> I, I)
            trimmed
                .chars()
                .map(|c| match c.to_ascii_uppercase() {
                    'A' | 'C' | 'G' | 'T' | 'N' => Allele::Base(c.to_string()),
                    'D' => Allele::Deletion,
                    'I' => Allele::Insertion,
                    '-' | '0' | '?' => Allele::Missing,
                    _ => Allele::Missing, // Treat garbage as missing
                })
                .collect()
        }
    }
}

/// Iterator over DTC records in a raw genotype text file.
pub struct Reader<R> {
    inner: R,
    line: u64,
    buf: String,
    has_warned_build: bool,
}

impl<R> Reader<R>
where
    R: BufRead,
{
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            line: 0,
            buf: String::new(),
            has_warned_build: false,
        }
    }

    pub fn into_inner(self) -> R {
        self.inner
    }
}

impl<R> Iterator for Reader<R>
where
    R: BufRead,
{
    type Item = Result<Record, ParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.buf.clear();
            match self.inner.read_line(&mut self.buf) {
                Ok(0) => return None,
                Ok(_) => {
                    self.line += 1;
                    // Sanitize input: remove quotes which are common in CSV formats (MyHeritage)
                    let trimmed = self.buf.trim_end_matches(&['\n', '\r'][..]).trim();

                    if trimmed.is_empty() {
                        continue;
                    }

                    if trimmed.starts_with('#') {
                        if !self.has_warned_build {
                            if trimmed.contains("Build 36") || trimmed.contains("NCBI36") {
                                tracing::warn!(
                                    "Input file appears to be Build 36/NCBI36! Coordinate mismatches with GRCh38 are likely."
                                );
                                self.has_warned_build = true;
                            } else if trimmed.contains("Build 37") || trimmed.contains("GRCh37") {
                                tracing::warn!(
                                    "Input file appears to be Build 37/GRCh37. Ensure you are using a compatible reference."
                                );
                                self.has_warned_build = true;
                            }
                        }
                        continue;
                    }

                    // Check for header by sanitizing quotes first (handling "RSID" in CSVs)
                    // Use safe slicing to avoid panic on multi-byte UTF-8 characters
                    let header_check = trimmed.trim_matches('"');
                    if let Some(prefix) = header_check.get(..4)
                        && (prefix.eq_ignore_ascii_case("rsid")
                            || prefix.eq_ignore_ascii_case("loid"))
                    {
                        continue;
                    }

                    return Some(parse_record(trimmed).map_err(|kind| ParseError {
                        line: self.line,
                        raw: trimmed.to_string(),
                        kind,
                    }));
                }
                Err(e) => {
                    return Some(Err(ParseError {
                        line: self.line,
                        raw: String::new(),
                        kind: ParseErrorKind::Io(e),
                    }));
                }
            }
        }
    }
}

/// Errors that can arise while parsing a DTC genotype record.
#[derive(Debug, Error)]
#[error("line {line}: {kind}")]
pub struct ParseError {
    pub line: u64,
    pub raw: String,
    #[source]
    pub kind: ParseErrorKind,
}

#[derive(Debug, Error)]
pub enum ParseErrorKind {
    #[error("I/O error")]
    Io(#[from] io::Error),
    #[error("expected 4, 5, or 6 fields, found {0}")]
    FieldCount(usize),
    #[error("invalid chromosome field")]
    InvalidChromosome,
    #[error("invalid position: {0}")]
    InvalidPosition(ParseIntError),
    #[error("missing genotype field")]
    MissingGenotype,
    #[error(
        "invalid variant ID {0:?}: must not contain whitespace or ';' (forbidden in the VCF ID field)"
    )]
    InvalidId(String),
}

/// Returns true if `id` is usable as a VCF ID field value.
///
/// Per the VCF spec (§1.6.1.3 "Fixed fields: ID"), the ID field must not
/// contain whitespace or semicolons. noodles enforces this at serialization
/// time and, when violated, fails deep in the external-sort spill writer with
/// a bare "invalid ID" error that does not name the offending line. We mirror
/// the check here so a malformed ID is rejected cleanly at parse time and the
/// offending record is skipped with a counted warning, consistent with the
/// other field-level parse errors.
fn is_valid_vcf_id(id: &str) -> bool {
    id.chars().all(|c| !c.is_whitespace() && c != ';')
}

fn parse_record(line: &str) -> Result<Record, ParseErrorKind> {
    // Determine delimiter (Comma for MyHeritage/FTDNA, Whitespace for others)
    // and strip quotes from fields lazily to avoid allocation
    let fields: Vec<&str> = if line.contains(',') {
        line.split(',')
            .map(|s| s.trim().trim_matches('"'))
            .collect()
    } else {
        line.split_whitespace()
            .map(|s| s.trim_matches('"'))
            .collect()
    };

    let count = fields.len();

    // 3. Heuristic column mapping based on field count
    let (id_idx, chr_idx, pos_idx, geno_val) = match count {
        4 => {
            // Standard / 23andMe
            // ID, Chr, Pos, Genotype
            (0, 1, 2, fields[3].to_string())
        }
        5 => {
            // AncestryDNA
            // ID, Chr, Pos, A1, A2
            (0, 1, 2, format!("{}/{}", fields[3], fields[4]))
        }
        6 => {
            // deCODEme
            // ID, Variation, Chr, Pos, Strand, Genotype
            let strand = fields[4];
            let mut g = fields[5].to_string();
            if strand == "-" {
                g = flip_genotype(&g);
            }
            (0, 2, 3, g)
        }
        _ => return Err(ParseErrorKind::FieldCount(count)),
    };

    let id_str = fields[id_idx];
    let chromosome = fields[chr_idx];
    let position_str = fields[pos_idx];

    if chromosome.is_empty() {
        return Err(ParseErrorKind::InvalidChromosome);
    }

    let position = position_str
        .parse::<u64>()
        .map_err(ParseErrorKind::InvalidPosition)?;

    let id = if id_str.is_empty() || id_str == "0" || id_str == "." {
        None
    } else if is_valid_vcf_id(id_str) {
        Some(id_str.to_string())
    } else {
        // The ID carries characters the VCF ID field forbids (whitespace or
        // ';'). If we kept it, noodles would later abort the whole conversion
        // with a bare "invalid ID" from inside the spill writer. Reject the
        // record here so it is skipped with a counted warning instead.
        return Err(ParseErrorKind::InvalidId(id_str.to_string()));
    };

    Ok(Record {
        id,
        chromosome: chromosome.to_string(),
        position,
        genotype: geno_val,
        metrics: None,
    })
}

fn flip_genotype(g: &str) -> String {
    g.chars()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'a' => 't',
            't' => 'a',
            'c' => 'g',
            'g' => 'c',
            other => other,
        })
        .collect()
}

impl fmt::Display for Record {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}",
            self.id.as_deref().unwrap_or("."),
            self.chromosome,
            self.position,
            self.genotype
        )
    }
}

/// Lossless on-disk spill codec for external sorting.
///
/// This is deliberately *separate* from the lenient input parser
/// ([`parse_record`]): the spill format is a fixed, self-describing 8-column
/// tab-delimited line that round-trips every field — including optional
/// [`ArrayMetrics`] — exactly. Decoupling it means the external sorter never
/// re-runs the field-count heuristics meant for heterogeneous DTC inputs, and
/// metrics survive a spill to disk.
///
/// Columns: `id, chr, pos, genotype, baf, lrr, gencall, gentrain`. A missing
/// id is `.`; a missing numeric metric is the empty string. When all four
/// metric columns are empty the record's `metrics` is reconstituted as `None`.
impl Record {
    pub fn to_spill_line(&self) -> String {
        fn num(v: Option<f32>) -> String {
            v.map(|x| x.to_string()).unwrap_or_default()
        }
        let m = self.metrics.unwrap_or_default();
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.id.as_deref().unwrap_or("."),
            self.chromosome,
            self.position,
            self.genotype,
            num(m.baf),
            num(m.lrr),
            num(m.gencall),
            num(m.gentrain),
        )
    }

    pub fn from_spill_line(line: &str) -> Option<Record> {
        // Strict: a spill line is exactly 8 tab-separated fields. Anything else
        // means the internal spill file is corrupt; reject rather than silently
        // reshaping the record (every field is tab-free by construction).
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() != 8 {
            return None;
        }
        let id_str = cols[0];
        let chromosome = cols[1];
        let position: u64 = cols[2].parse().ok()?;
        let genotype = cols[3];

        fn num(s: &str) -> Option<f32> {
            if s.is_empty() {
                None
            } else {
                s.parse::<f32>().ok().filter(|v| v.is_finite())
            }
        }
        let metrics = ArrayMetrics {
            baf: num(cols[4]),
            lrr: num(cols[5]),
            gencall: num(cols[6]),
            gentrain: num(cols[7]),
        };

        Some(Record {
            id: if id_str == "." || id_str.is_empty() {
                None
            } else {
                Some(id_str.to_string())
            },
            chromosome: chromosome.to_string(),
            position,
            genotype: genotype.to_string(),
            metrics: if metrics.is_empty() {
                None
            } else {
                Some(metrics)
            },
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_basic_record() {
        let record = parse_record("rs1\t1\t42\tAG").expect("parse");
        assert_eq!(record.id.as_deref(), Some("rs1"));
        assert_eq!(record.chromosome, "1");
        assert_eq!(record.position, 42);
        assert_eq!(record.genotype, "AG");
        assert_eq!(record.metrics, None);
    }

    #[test]
    fn spill_line_round_trips_without_metrics() {
        let rec = Record {
            id: Some("rs1".into()),
            chromosome: "X".into(),
            position: 12345,
            genotype: "AG".into(),
            metrics: None,
        };
        let back = Record::from_spill_line(&rec.to_spill_line()).expect("round-trip");
        assert_eq!(back, rec);
        // A missing id round-trips to None, too.
        let anon = Record {
            id: None,
            chromosome: "1".into(),
            position: 1,
            genotype: "--".into(),
            metrics: None,
        };
        assert_eq!(Record::from_spill_line(&anon.to_spill_line()), Some(anon));
    }

    #[test]
    fn spill_line_round_trips_metrics() {
        let rec = Record {
            id: Some("rs9".into()),
            chromosome: "12".into(),
            position: 88675,
            genotype: "CT".into(),
            metrics: Some(ArrayMetrics {
                baf: Some(0.4653),
                lrr: Some(-0.0240),
                gencall: Some(0.8985),
                gentrain: Some(0.8622),
            }),
        };
        let back = Record::from_spill_line(&rec.to_spill_line()).expect("round-trip");
        assert_eq!(back, rec);
    }

    #[test]
    fn spill_line_round_trips_partial_metrics() {
        let rec = Record {
            id: Some("rs9".into()),
            chromosome: "12".into(),
            position: 88675,
            genotype: "CT".into(),
            metrics: Some(ArrayMetrics {
                baf: Some(0.5),
                lrr: None,
                gencall: None,
                gentrain: None,
            }),
        };
        let back = Record::from_spill_line(&rec.to_spill_line()).expect("round-trip");
        assert_eq!(back, rec);
    }

    #[test]
    fn reader_skips_comments_and_detects_build() {
        let data = b"#comment\n# Build 37\nrs1\t1\t10\tAA\n";
        let mut reader = Reader::new(&data[..]);
        let record = reader.next().unwrap().unwrap();
        assert_eq!(record.position, 10);
        assert!(reader.has_warned_build);
        assert!(reader.next().is_none());
    }

    #[test]
    fn parses_csv_quoted() {
        // MyHeritage style
        let line = "\"rs123\",\"1\",\"100\",\"AA\"";
        let record = parse_record(line).unwrap();
        assert_eq!(record.id.as_deref(), Some("rs123"));
        assert_eq!(record.chromosome, "1");
        assert_eq!(record.position, 100);
        assert_eq!(record.genotype, "AA");
    }

    #[test]
    fn parses_ancestry_5col() {
        // AncestryDNA style: ID, Chr, Pos, A1, A2
        let line = "rs123\t1\t100\tA\tG";
        let record = parse_record(line).unwrap();
        assert_eq!(record.genotype, "A/G");
    }

    #[test]
    fn rejects_id_with_semicolon() {
        // An ID carrying a ';' is forbidden in the VCF ID field. Previously
        // this survived parsing and crashed the whole conversion deep in the
        // external-sort spill writer with a bare "invalid ID". It must now be
        // rejected cleanly at parse time so the record is skipped with a
        // counted warning (like the other field-level parse errors).
        let err = parse_record("rs1;weird\t1\t100\tAA").unwrap_err();
        match &err {
            ParseErrorKind::InvalidId(id) => assert_eq!(id, "rs1;weird"),
            other => panic!("expected InvalidId, got {other:?}"),
        }
    }

    #[test]
    fn rejects_id_with_embedded_whitespace() {
        // Tab/newline-split fields can't carry spaces, but other Unicode
        // whitespace (e.g. a no-break space) can sneak through a CSV field and
        // is equally rejected by the VCF ID grammar.
        let err = parse_record("rs1\u{00a0}x,1,100,AA").unwrap_err();
        assert!(matches!(err, ParseErrorKind::InvalidId(_)));
    }

    #[test]
    fn accepts_normal_rsid() {
        // Regression guard: ordinary IDs (and the missing-ID sentinels) must
        // still parse fine after adding ID validation.
        assert_eq!(
            parse_record("rs123\t1\t100\tAA").unwrap().id.as_deref(),
            Some("rs123")
        );
        assert_eq!(parse_record(".\t1\t100\tAA").unwrap().id, None);
        // Colons and other punctuation that the VCF ID grammar allows are kept
        // (synthetic IDs use ':' as a separator).
        assert_eq!(
            parse_record("1:100:A:G\t1\t100\tAA")
                .unwrap()
                .id
                .as_deref(),
            Some("1:100:A:G")
        );
    }

    #[test]
    fn parses_decodeme_6col_flipped() {
        // deCODEme style: ID, Var, Chr, Pos, Strand, Geno
        // Strand - means flip A->T
        let line = "rs123\tvar\t1\t100\t-\tA";
        let record = parse_record(line).unwrap();
        assert_eq!(record.genotype, "T"); // Flipped from A
    }
}
