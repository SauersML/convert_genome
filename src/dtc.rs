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
}

/// Iterator over DTC records in a raw genotype text file.
pub struct Reader<R> {
    inner: R,
    line: u64,
    buf: String,
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
                    let trimmed = self.buf.trim_end_matches(&['\n', '\r'][..]);
                    if trimmed.is_empty() || trimmed.starts_with('#') {
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
    #[error("expected at least four tab- or space-delimited fields, found {0}")]
    FieldCount(usize),
    #[error("invalid chromosome field")]
    InvalidChromosome,
    #[error("invalid position: {0}")]
    InvalidPosition(ParseIntError),
    #[error("missing genotype field")]
    MissingGenotype,
}

fn parse_record(line: &str) -> Result<Record, ParseErrorKind> {
    let mut fields = line.split_whitespace();

    let id_field = fields.next();
    let chromosome = fields
        .next()
        .ok_or_else(|| ParseErrorKind::FieldCount(count_fields(line)))?;
    if chromosome.is_empty() {
        return Err(ParseErrorKind::InvalidChromosome);
    }

    let position_str = fields
        .next()
        .ok_or_else(|| ParseErrorKind::FieldCount(count_fields(line)))?;
    let position = position_str
        .parse::<u64>()
        .map_err(ParseErrorKind::InvalidPosition)?;

    let genotype = fields
        .next()
        .ok_or(ParseErrorKind::MissingGenotype)?
        .to_string();

    let id = id_field.and_then(|s| {
        let trimmed = s.trim();
        if trimmed.is_empty() || trimmed == "0" || trimmed == "." {
            None
        } else {
            Some(trimmed.to_string())
        }
    });

    Ok(Record {
        id,
        chromosome: chromosome.to_string(),
        position,
        genotype,
    })
}

fn count_fields(line: &str) -> usize {
    line.split_whitespace().count()
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
    fn reader_skips_comments() {
        let data = b"#comment\nrs1\t1\t10\tAA\n";
        let mut reader = Reader::new(&data[..]);
        let record = reader.next().unwrap().unwrap();
        assert_eq!(record.position, 10);
        assert!(reader.next().is_none());
    }
}
