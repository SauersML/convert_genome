//! Reader for Illumina GenomeStudio Genotyping ("GSGT") Final Report exports.
//!
//! A Final Report is the standard genotype+intensity export from an Illumina
//! Infinium array (GSA, OmniExpress, etc.). Its layout is:
//!
//! ```text
//! [Header]
//! GSGT Version,2.0.5
//! Processing Date,5/19/2025 12:53 PM
//! Content,,SomeManifest.bpm
//! Num SNPs,730528
//! Num Samples,1
//! [Data]
//! Sample ID,SNP Name,Chr,Position,...,Allele1 - Plus,Allele2 - Plus,...
//! GRC123,rs123,1,102914837,...,G,G,...
//! ...
//! ```
//!
//! Unlike the simple DTC text formats (23andMe etc.), the column set is
//! *configurable* in GenomeStudio, so we map columns by name rather than
//! by position.
//!
//! ## Strand
//!
//! GenomeStudio can emit alleles in four conventions: `Plus`, `Top`,
//! `Forward`, and `AB`. Only `Plus` is consistently on the genome's `+`
//! strand for every marker — `Top`/`Forward`/`AB` flip per-SNP. Because the
//! conversion pipeline models strand *file-wide*, we deliberately require the
//! `Plus` alleles: they are authoritative forward-strand calls, so no strand
//! inference is needed and no per-SNP corruption can occur. A report exported
//! without Plus alleles is rejected with a clear error rather than guessed.

use std::io::{self, BufRead};

use thiserror::Error;

use crate::dtc::{ArrayMetrics, Record as DtcRecord};

/// Metadata captured from the `[Header]` block of a Final Report.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ReportMetadata {
    /// `GSGT Version` value, e.g. `2.0.5`.
    pub gsgt_version: Option<String>,
    /// The manifest (`.bpm`) the chip was processed against, from `Content`.
    pub manifest: Option<String>,
    /// `Num SNPs` value, if present and parseable.
    pub num_snps: Option<u64>,
    /// `Num Samples` value, if present and parseable.
    pub num_samples: Option<u64>,
}

/// Errors that can arise while reading a GenomeStudio Final Report.
#[derive(Debug, Error)]
pub enum GsError {
    #[error("I/O error")]
    Io(#[from] io::Error),
    #[error("not a GenomeStudio Final Report: no [Data] section found")]
    MissingDataSection,
    #[error("GenomeStudio report has no column header row after [Data]")]
    MissingColumnHeader,
    #[error("GenomeStudio report is missing required column(s): {0}")]
    MissingColumns(String),
    #[error(
        "GenomeStudio report has no 'Allele1/2 - Plus' columns. Re-export the Final \
         Report with Plus-strand alleles enabled; Top/Forward/AB cannot be converted \
         accurately without the array manifest."
    )]
    NoPlusAlleles,
    #[error("line {line}: invalid position '{value}'")]
    InvalidPosition { line: u64, value: String },
    #[error(
        "line {line}: invalid variant ID {value:?}: must not contain whitespace or ';' \
         (forbidden in the VCF ID field)"
    )]
    InvalidId { line: u64, value: String },
    #[error("line {line}: wrong field count (expected {expected}, found {found})")]
    FieldCount {
        line: u64,
        expected: usize,
        found: usize,
    },
}

/// Resolved column indices for the fields we consume.
#[derive(Debug, Clone, Copy)]
struct ColumnMap {
    sample_id: Option<usize>,
    snp_name: usize,
    chr: usize,
    position: usize,
    allele1_plus: usize,
    allele2_plus: usize,
    /// Optional per-sample metric columns, carried into VCF FORMAT fields.
    baf: Option<usize>,
    lrr: Option<usize>,
    gencall: Option<usize>,
    gentrain: Option<usize>,
    /// Number of columns in the header; used to validate row width.
    width: usize,
}

/// Normalize a column name for matching: trim, strip surrounding quotes, and
/// lowercase. We keep internal spacing so `"allele1 - plus"` matches.
fn norm(name: &str) -> String {
    name.trim().trim_matches('"').trim().to_ascii_lowercase()
}

impl ColumnMap {
    fn from_header(header_fields: &[&str]) -> Result<Self, GsError> {
        let normalized: Vec<String> = header_fields.iter().map(|f| norm(f)).collect();

        let find = |candidates: &[&str]| -> Option<usize> {
            candidates
                .iter()
                .find_map(|c| normalized.iter().position(|n| n == c))
        };
        // Looser match: any column whose normalized name ends with `suffix`.
        let find_suffix =
            |suffix: &str| -> Option<usize> { normalized.iter().position(|n| n.ends_with(suffix)) };

        // SNP name: prefer the canonical names, then anything ending in "snp name"
        // (e.g. the "GxG SNP Name" variant seen on some GSA exports).
        let snp_name = find(&["snp name", "gxg snp name", "name"])
            .or_else(|| find_suffix("snp name"))
            .ok_or_else(|| GsError::MissingColumns("SNP Name".into()))?;

        let chr =
            find(&["chr", "chromosome"]).ok_or_else(|| GsError::MissingColumns("Chr".into()))?;

        let position = find(&["position", "pos", "mapinfo"])
            .ok_or_else(|| GsError::MissingColumns("Position".into()))?;

        let allele1_plus = find(&["allele1 - plus"]);
        let allele2_plus = find(&["allele2 - plus"]);
        let (allele1_plus, allele2_plus) = match (allele1_plus, allele2_plus) {
            (Some(a), Some(b)) => (a, b),
            _ => return Err(GsError::NoPlusAlleles),
        };

        Ok(ColumnMap {
            sample_id: find(&["sample id", "sample_id", "sample"]),
            snp_name,
            chr,
            position,
            allele1_plus,
            allele2_plus,
            baf: find(&["b allele freq", "baf", "b allele frequency"]),
            lrr: find(&["log r ratio", "lrr"]),
            gencall: find(&["gc score", "gencall score", "gencall"]),
            gentrain: find(&["gt score", "gentrain score", "gentrain"]),
            width: header_fields.len(),
        })
    }
}

/// Streaming reader over the `[Data]` rows of a GenomeStudio Final Report.
///
/// Yields [`DtcRecord`]s so the existing conversion pipeline can consume the
/// report with no special-casing. The genotype is built from the two
/// `Plus`-strand alleles.
pub struct Reader<R> {
    inner: R,
    columns: ColumnMap,
    line: u64,
    buf: String,
    /// The first `Sample ID` we observed; rows from other samples are skipped.
    locked_sample: Option<String>,
    /// Count of rows skipped because they belonged to a second/third sample.
    pub skipped_other_sample: u64,
    metadata: ReportMetadata,
}

impl<R> Reader<R>
where
    R: BufRead,
{
    /// Build a reader, consuming the `[Header]` block and column-header row.
    ///
    /// Returns an error if the stream is not a Final Report or lacks the
    /// required columns (including `Plus`-strand alleles).
    pub fn new(mut inner: R) -> Result<Self, GsError> {
        let mut metadata = ReportMetadata::default();
        let mut line: u64 = 0;
        let mut buf = String::new();

        // 1. Walk the [Header] block, capturing metadata, until [Data].
        loop {
            buf.clear();
            let n = inner.read_line(&mut buf)?;
            if n == 0 {
                return Err(GsError::MissingDataSection);
            }
            line += 1;
            let trimmed = strip_bom(buf.trim_end_matches(['\n', '\r']).trim());

            if trimmed.is_empty() {
                continue;
            }
            if trimmed.eq_ignore_ascii_case("[data]") {
                break;
            }
            if trimmed.eq_ignore_ascii_case("[header]") {
                continue;
            }
            // Header metadata lines are "Key,Value" (Content uses "Key,,Value").
            absorb_metadata(trimmed, &mut metadata);
        }

        // 2. Next non-empty line is the column header.
        let header_line = loop {
            buf.clear();
            let n = inner.read_line(&mut buf)?;
            if n == 0 {
                return Err(GsError::MissingColumnHeader);
            }
            line += 1;
            let trimmed = buf.trim_end_matches(['\n', '\r']).trim();
            if !trimmed.is_empty() {
                break trimmed.to_string();
            }
        };

        let header_fields: Vec<&str> = header_line.split(',').collect();
        let columns = ColumnMap::from_header(&header_fields)?;

        Ok(Self {
            inner,
            columns,
            line,
            buf: String::new(),
            locked_sample: None,
            skipped_other_sample: 0,
            metadata,
        })
    }

    /// Metadata parsed from the `[Header]` block.
    pub fn metadata(&self) -> &ReportMetadata {
        &self.metadata
    }

    fn parse_data_line(&mut self, raw: &str) -> Result<Option<DtcRecord>, GsError> {
        let fields: Vec<&str> = raw.split(',').map(|s| s.trim()).collect();
        // Require the exact column count. A short row means truncation; a long
        // row means an embedded delimiter shifted every column (GenomeStudio
        // Final Reports are unquoted CSV, so a stray comma silently misaligns
        // SNP/chr/pos/genotype). Either way we refuse to emit a misread record.
        // A single trailing empty field (stray terminal comma) is tolerated.
        let n = fields.len();
        let width = self.columns.width;
        let acceptable = n == width || (n == width + 1 && fields[n - 1].is_empty());
        if !acceptable {
            return Err(GsError::FieldCount {
                line: self.line,
                expected: width,
                found: n,
            });
        }

        // Single-sample lock: convert_genome emits one sample. Process the
        // first Sample ID we see; skip rows belonging to later samples.
        if let Some(idx) = self.columns.sample_id {
            let sample = fields[idx];
            match &self.locked_sample {
                None => self.locked_sample = Some(sample.to_string()),
                Some(locked) if locked != sample => {
                    self.skipped_other_sample += 1;
                    return Ok(None);
                }
                Some(_) => {}
            }
        }

        let chromosome = fields[self.columns.chr].trim();
        if chromosome.is_empty() {
            // No useful coordinate; skip silently (intensity-only probes).
            return Ok(None);
        }

        let position_str = fields[self.columns.position].trim();
        let position = position_str
            .parse::<u64>()
            .map_err(|_| GsError::InvalidPosition {
                line: self.line,
                value: position_str.to_string(),
            })?;
        if position == 0 {
            // Unmapped probe (Position 0); nothing to emit.
            return Ok(None);
        }

        let a1 = fields[self.columns.allele1_plus].trim();
        let a2 = fields[self.columns.allele2_plus].trim();
        let genotype = join_alleles(a1, a2);

        let snp_name = fields[self.columns.snp_name].trim();
        let id = if snp_name.is_empty() || snp_name == "." {
            None
        } else if crate::dtc::is_valid_vcf_id(snp_name) {
            Some(snp_name.to_string())
        } else {
            // A SNP name carrying whitespace or ';' is forbidden in the VCF ID
            // field. Left intact it would survive into the RecordBuf and later
            // abort the whole conversion deep in the serializer with a bare
            // "invalid ID". Reject it here so the row is skipped with a counted
            // warning, mirroring the DTC parser.
            return Err(GsError::InvalidId {
                line: self.line,
                value: snp_name.to_string(),
            });
        };

        let num = |idx: Option<usize>| -> Option<f32> {
            let s = fields[idx?].trim();
            // Empty and non-finite values (GenomeStudio writes `NaN` for
            // intensity-only loci) become missing, not invalid VCF floats.
            if s.is_empty() {
                return None;
            }
            s.parse::<f32>().ok().filter(|v| v.is_finite())
        };
        let metrics = ArrayMetrics {
            baf: num(self.columns.baf),
            lrr: num(self.columns.lrr),
            gencall: num(self.columns.gencall),
            gentrain: num(self.columns.gentrain),
        };

        Ok(Some(DtcRecord {
            id,
            chromosome: chromosome.to_string(),
            position,
            genotype,
            metrics: if metrics.is_empty() {
                None
            } else {
                Some(metrics)
            },
        }))
    }
}

impl<R> Iterator for Reader<R>
where
    R: BufRead,
{
    type Item = Result<DtcRecord, GsError>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.buf.clear();
            match self.inner.read_line(&mut self.buf) {
                Ok(0) => return None,
                Ok(_) => {
                    self.line += 1;
                    let trimmed = self.buf.trim_end_matches(['\n', '\r']).trim();
                    if trimmed.is_empty() {
                        continue;
                    }
                    // Some exports append a trailing "[Data]"-less footer? None
                    // known, but guard against stray section markers.
                    if trimmed.starts_with('[') && trimmed.ends_with(']') {
                        continue;
                    }
                    // Take ownership so we don't borrow self.buf while &mut self.
                    let owned = trimmed.to_string();
                    match self.parse_data_line(&owned) {
                        Ok(Some(rec)) => return Some(Ok(rec)),
                        Ok(None) => continue,
                        Err(e) => return Some(Err(e)),
                    }
                }
                Err(e) => return Some(Err(GsError::Io(e))),
            }
        }
    }
}

/// Join the two Plus-strand alleles into a genotype string the DTC allele
/// parser understands. Single-character alleles concatenate (`"AG"`, `"--"`);
/// multi-character alleles (insertion sequences) are `/`-separated.
fn join_alleles(a1: &str, a2: &str) -> String {
    let a1 = if a1.is_empty() { "-" } else { a1 };
    let a2 = if a2.is_empty() { "-" } else { a2 };
    if a1.chars().count() <= 1 && a2.chars().count() <= 1 {
        format!("{a1}{a2}")
    } else {
        format!("{a1}/{a2}")
    }
}

/// Parse a `"Key,Value"` (or `"Key,,Value"`) header line into metadata.
fn absorb_metadata(line: &str, meta: &mut ReportMetadata) {
    let mut parts = line.splitn(2, ',');
    let key = parts.next().unwrap_or("").trim();
    let rest = parts.next().unwrap_or("").trim();
    // `Content` uses an empty middle field: "Content,,Manifest.bpm".
    let value = rest.trim_start_matches(',').trim();
    match key.to_ascii_lowercase().as_str() {
        "gsgt version" => meta.gsgt_version = Some(value.to_string()),
        "content" => meta.manifest = Some(value.to_string()),
        "num snps" => meta.num_snps = value.parse().ok(),
        "num samples" => meta.num_samples = value.parse().ok(),
        _ => {}
    }
}

fn strip_bom(s: &str) -> &str {
    s.strip_prefix('\u{feff}').unwrap_or(s)
}

/// Heuristic detector: does this leading chunk look like a GenomeStudio Final
/// Report? True when the first non-empty, non-BOM content is the `[Header]`
/// (or `[Data]`) section marker, or mentions `GSGT Version`.
pub fn looks_like_genome_studio(buf: &[u8]) -> bool {
    // Inspect only the first few non-empty lines.
    let text = match std::str::from_utf8(buf) {
        Ok(t) => t,
        Err(_) => {
            // Fall back to lossy on the valid prefix.
            let valid = buf.iter().take_while(|&&b| b.is_ascii()).count();
            std::str::from_utf8(&buf[..valid]).unwrap_or("")
        }
    };
    for line in text.lines().take(8) {
        let t = strip_bom(line.trim());
        if t.is_empty() {
            continue;
        }
        if t.eq_ignore_ascii_case("[header]") || t.eq_ignore_ascii_case("[data]") {
            return true;
        }
        if t.to_ascii_lowercase().starts_with("gsgt version") {
            return true;
        }
        // First meaningful line is something else → not a Final Report.
        return false;
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "\
[Header]
GSGT Version,2.0.5
Processing Date,5/19/2025 12:53 PM
Content,,GxGComprehensiveGSAv2-2_20031858_A2.bpm
Num SNPs,4
Num Samples,1
[Data]
Sample ID,SNP Name,Chr,Position,Log R Ratio,B Allele Freq,Allele1 - Plus,Allele2 - Plus,Allele1 - Top,Allele2 - Top
S1,rs100,1,102914837,-0.02,1.00,G,G,G,G
S1,rs200,1,106194696,-0.24,0.07,A,T,A,A
S1,rs300,2,500,0.00,0.50,-,-,A,G
S1,rs400,2,0,0.00,0.50,C,C,C,C
";

    fn reader(s: &str) -> Reader<&[u8]> {
        Reader::new(s.as_bytes()).expect("construct reader")
    }

    #[test]
    fn parses_header_metadata() {
        let r = reader(SAMPLE);
        let m = r.metadata();
        assert_eq!(m.gsgt_version.as_deref(), Some("2.0.5"));
        assert_eq!(
            m.manifest.as_deref(),
            Some("GxGComprehensiveGSAv2-2_20031858_A2.bpm")
        );
        assert_eq!(m.num_snps, Some(4));
        assert_eq!(m.num_samples, Some(1));
    }

    #[test]
    fn extracts_plus_alleles_and_skips_unmapped() {
        let recs: Vec<_> = reader(SAMPLE).map(|r| r.unwrap()).collect();
        // rs400 has Position 0 (unmapped) → skipped.
        assert_eq!(recs.len(), 3);

        assert_eq!(recs[0].id.as_deref(), Some("rs100"));
        assert_eq!(recs[0].chromosome, "1");
        assert_eq!(recs[0].position, 102914837);
        assert_eq!(recs[0].genotype, "GG");

        // Heterozygous, taken from Plus (A/T), not Top (A/A).
        assert_eq!(recs[1].genotype, "AT");

        // No-call "-"/"-" → "--", which DtcRecord::is_missing recognizes.
        assert_eq!(recs[2].genotype, "--");
        assert!(recs[2].is_missing());
    }

    #[test]
    fn parses_array_metrics() {
        // recs[0] = rs100 with Log R Ratio -0.02, B Allele Freq 1.00; the
        // SAMPLE report has no GC/GT Score columns, so those stay None.
        let recs: Vec<_> = reader(SAMPLE).map(|r| r.unwrap()).collect();
        let m = recs[0].metrics.expect("metrics present");
        assert_eq!(m.lrr, Some(-0.02));
        assert_eq!(m.baf, Some(1.00));
        assert_eq!(m.gencall, None);
        assert_eq!(m.gentrain, None);

        // A report carrying all four metric columns populates all four.
        let full = "\
[Data]
SNP Name,Chr,Position,Log R Ratio,B Allele Freq,GC Score,GT Score,Allele1 - Plus,Allele2 - Plus
rs1,1,100,-0.12,0.48,0.97,0.81,A,G
";
        let r = Reader::new(full.as_bytes())
            .unwrap()
            .next()
            .unwrap()
            .unwrap();
        let m = r.metrics.expect("metrics present");
        assert_eq!(m.lrr, Some(-0.12));
        assert_eq!(m.baf, Some(0.48));
        assert_eq!(m.gencall, Some(0.97));
        assert_eq!(m.gentrain, Some(0.81));
    }

    #[test]
    fn rejects_snp_name_with_forbidden_chars() {
        // A SNP name carrying a ';' (or whitespace) is forbidden in the VCF ID
        // field. It must be rejected at parse time as GsError::InvalidId so the
        // row is skipped with a counted warning, instead of surviving into the
        // RecordBuf and crashing the whole conversion deep in the serializer.
        let report = "\
[Data]
SNP Name,Chr,Position,Allele1 - Plus,Allele2 - Plus
rs1;weird,1,100,A,G
rs2,1,200,C,T
";
        let results: Vec<_> = Reader::new(report.as_bytes()).unwrap().collect();
        assert_eq!(results.len(), 2);
        match &results[0] {
            Err(GsError::InvalidId { value, .. }) => assert_eq!(value, "rs1;weird"),
            other => panic!("expected InvalidId, got {other:?}"),
        }
        // The well-formed row still parses fine.
        assert_eq!(results[1].as_ref().unwrap().id.as_deref(), Some("rs2"));
    }

    #[test]
    fn requires_plus_alleles() {
        let no_plus = "\
[Header]
GSGT Version,2.0.5
[Data]
SNP Name,Chr,Position,Allele1 - Top,Allele2 - Top
rs1,1,100,A,G
";
        let err = Reader::new(no_plus.as_bytes()).err();
        assert!(
            matches!(err, Some(GsError::NoPlusAlleles)),
            "expected NoPlusAlleles, got {err:?}"
        );
    }

    #[test]
    fn locks_onto_first_sample() {
        let multi = "\
[Header]
Num Samples,2
[Data]
Sample ID,SNP Name,Chr,Position,Allele1 - Plus,Allele2 - Plus
S1,rs1,1,100,A,A
S2,rs1,1,100,C,C
S1,rs2,1,200,G,T
S2,rs2,1,200,T,T
";
        let mut r = Reader::new(multi.as_bytes()).unwrap();
        let recs: Vec<_> = (&mut r).map(|x| x.unwrap()).collect();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].genotype, "AA");
        assert_eq!(recs[1].genotype, "GT");
        assert_eq!(r.skipped_other_sample, 2);
    }

    #[test]
    fn column_order_is_by_name_not_position() {
        // Columns deliberately reordered; GxG SNP Name variant; extra columns.
        let reordered = "\
[Data]
Chr,Position,Allele2 - Plus,GxG SNP Name,GC Score,Allele1 - Plus,Sample ID
1,100,T,rsX,0.99,A,S1
";
        let recs: Vec<_> = Reader::new(reordered.as_bytes())
            .unwrap()
            .map(|x| x.unwrap())
            .collect();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].id.as_deref(), Some("rsX"));
        assert_eq!(recs[0].chromosome, "1");
        assert_eq!(recs[0].position, 100);
        // Allele1 - Plus = A, Allele2 - Plus = T → "AT".
        assert_eq!(recs[0].genotype, "AT");
    }

    #[test]
    fn detector_recognizes_reports() {
        assert!(looks_like_genome_studio(b"[Header]\nGSGT Version,2.0.5\n"));
        assert!(looks_like_genome_studio(b"\n\n[Data]\nSNP Name,Chr\n"));
        assert!(looks_like_genome_studio(
            b"GSGT Version,2.0.5\nProcessing Date,x\n"
        ));
        // 23andMe / DTC files must NOT be detected as GenomeStudio.
        assert!(!looks_like_genome_studio(b"# rsid\trs1\t1\t100\tAA\n"));
        assert!(!looks_like_genome_studio(b"rsid\tchromosome\tposition\n"));
    }
}
