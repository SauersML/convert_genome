use std::{
    fmt,
    io::{self, BufRead},
    num::ParseIntError,
};

use thiserror::Error;

/// A single genotype entry from a direct-to-consumer text export.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Record {
    pub id: Option<String>,
    pub chromosome: String,
    pub position: u64,
    pub genotype: String,
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
            .map(|s| match s.trim() {
                "D" => Allele::Deletion,
                "I" => Allele::Insertion,
                "-" | "0" | "?" => Allele::Missing,
                val => Allele::Base(val.to_string()),
            })
            .collect()
    } else {
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
                                tracing::warn!("Input file appears to be Build 36/NCBI36! Coordinate mismatches with GRCh38 are likely.");
                                self.has_warned_build = true;
                            } else if trimmed.contains("Build 37") || trimmed.contains("GRCh37") {
                                tracing::warn!("Input file appears to be Build 37/GRCh37. Ensure you are using a compatible reference.");
                                self.has_warned_build = true;
                            }
                        }
                        continue;
                    }

                    // Check for header by sanitizing quotes first (handling "RSID" in CSVs)
                    let header_check = trimmed.trim_matches('"');
                    if header_check.len() >= 4 && (header_check[..4].eq_ignore_ascii_case("rsid") || header_check[..4].eq_ignore_ascii_case("loid")) {
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
}

fn parse_record(line: &str) -> Result<Record, ParseErrorKind> {
    // Determine delimiter (Comma for MyHeritage/FTDNA, Whitespace for others)
    // and strip quotes from fields lazily to avoid allocation
    let fields: Vec<&str> = if line.contains(',') {
        line.split(',').map(|s| s.trim().trim_matches('"')).collect()
    } else {
        line.split_whitespace().map(|s| s.trim_matches('"')).collect()
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
    } else {
        Some(id_str.to_string())
    };

    Ok(Record {
        id,
        chromosome: chromosome.to_string(),
        position,
        genotype: geno_val,
    })
}

fn flip_genotype(g: &str) -> String {
    g.chars()
        .map(|c| match c {
            'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C',
            'a' => 't', 't' => 'a', 'c' => 'g', 'g' => 'c',
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
    fn parses_decodeme_6col_flipped() {
        // deCODEme style: ID, Var, Chr, Pos, Strand, Geno
        // Strand - means flip A->T
        let line = "rs123\tvar\t1\t100\t-\tA";
        let record = parse_record(line).unwrap();
        assert_eq!(record.genotype, "T"); // Flipped from A
    }    
}
