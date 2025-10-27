use std::io::{BufRead, BufReader, Read};

use anyhow::{Context, Result, anyhow};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DtcRecord {
    pub id: String,
    pub chrom: String,
    pub position: u64,
    pub genotype: DtcGenotype,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DtcGenotype {
    pub alleles: [Allele; 2],
    pub is_call: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Allele {
    A,
    C,
    G,
    T,
    Deletion,
    Insertion,
    Other(char),
    Missing,
}

impl Allele {
    pub fn from_char(c: char) -> Self {
        match c {
            'A' | 'a' => Self::A,
            'C' | 'c' => Self::C,
            'G' | 'g' => Self::G,
            'T' | 't' => Self::T,
            '-' => Self::Missing,
            'D' | 'd' => Self::Deletion,
            'I' | 'i' => Self::Insertion,
            other => Self::Other(other),
        }
    }

    pub fn to_base_char(self) -> Option<char> {
        match self {
            Allele::A => Some('A'),
            Allele::C => Some('C'),
            Allele::G => Some('G'),
            Allele::T => Some('T'),
            Allele::Deletion | Allele::Insertion | Allele::Other(_) | Allele::Missing => None,
        }
    }
}

impl DtcGenotype {
    pub fn missing() -> Self {
        Self {
            alleles: [Allele::Missing; 2],
            is_call: false,
        }
    }

    pub fn from_call(call: &str) -> Self {
        if call.eq("--") || call.trim().is_empty() {
            return Self::missing();
        }

        let mut chars = call.chars().filter(|c| !c.is_whitespace());
        let first = chars
            .next()
            .map(Allele::from_char)
            .unwrap_or(Allele::Missing);
        let second = chars.next().map(Allele::from_char).unwrap_or(first);
        let is_call = !(matches!(first, Allele::Missing) && matches!(second, Allele::Missing));
        Self {
            alleles: [first, second],
            is_call,
        }
    }

    pub fn allele_chars(&self) -> [Option<char>; 2] {
        [
            self.alleles[0].to_base_char(),
            self.alleles[1].to_base_char(),
        ]
    }

    pub fn is_missing(&self) -> bool {
        !self.is_call
    }
}

pub struct DtcReader<R: Read> {
    inner: BufReader<R>,
}

impl<R: Read> DtcReader<R> {
    pub fn new(inner: R) -> Self {
        Self {
            inner: BufReader::new(inner),
        }
    }

    pub fn records(mut self) -> impl Iterator<Item = Result<DtcRecord>> {
        std::iter::from_fn(move || {
            loop {
                let mut line = String::new();
                match self.inner.read_line(&mut line) {
                    Ok(0) => return None,
                    Ok(_) => {
                        if line.trim().is_empty() || line.starts_with('#') {
                            continue;
                        }

                        match parse_record_line(line.trim_end()) {
                            Ok(Some(record)) => return Some(Ok(record)),
                            Ok(None) => continue,
                            Err(e) => return Some(Err(e)),
                        }
                    }
                    Err(e) => return Some(Err(anyhow!(e))),
                }
            }
        })
    }
}

fn normalize_chromosome(raw: &str) -> Option<String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return None;
    }

    let normalized = trimmed.trim_start_matches("chr");
    let normalized = normalized.trim_start_matches("CHR");
    let normalized = normalized.trim();

    let result = match normalized {
        "M" | "MT" | "m" | "mt" => "chrM".to_string(),
        "XY" | "xy" => "chrXY".to_string(),
        value if value.eq_ignore_ascii_case("PAR1") => "chrPAR1".to_string(),
        value if value.eq_ignore_ascii_case("PAR2") => "chrPAR2".to_string(),
        value => {
            if value.chars().all(|c| c.is_ascii_digit()) {
                match value.parse::<u32>() {
                    Ok(23) => "chrX".to_string(),
                    Ok(24) => "chrY".to_string(),
                    Ok(25) => "chrM".to_string(),
                    Ok(n) => format!("chr{}", n),
                    Err(_) => format!("chr{}", value),
                }
            } else if value.len() == 1
                && matches!(value.chars().next().unwrap(), 'X' | 'Y' | 'x' | 'y')
            {
                format!("chr{}", value.to_uppercase())
            } else {
                format!("chr{}", value)
            }
        }
    };

    Some(result)
}

fn parse_record_line(line: &str) -> Result<Option<DtcRecord>> {
    if line.is_empty() || line.starts_with('#') {
        return Ok(None);
    }

    let mut fields = line.split('\t');
    let id = fields
        .next()
        .ok_or_else(|| anyhow!("missing id column"))?
        .trim()
        .to_owned();
    let chrom = fields
        .next()
        .ok_or_else(|| anyhow!("missing chromosome column"))?
        .trim()
        .to_owned();
    let position = fields
        .next()
        .ok_or_else(|| anyhow!("missing position column"))?
        .trim()
        .parse::<u64>()
        .with_context(|| format!("invalid position for {id}"))?;
    let genotype = fields
        .next()
        .ok_or_else(|| anyhow!("missing genotype column"))?
        .trim()
        .to_owned();

    if fields.next().is_some() {
        return Err(anyhow!("unexpected additional columns in DTC record"));
    }

    let chrom = normalize_chromosome(&chrom).ok_or_else(|| anyhow!("invalid chromosome value"))?;
    let genotype = DtcGenotype::from_call(&genotype);

    Ok(Some(DtcRecord {
        id,
        chrom,
        position,
        genotype,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_basic_record() {
        let record = parse_record_line("rs1\t1\t101\tAG").unwrap().unwrap();
        assert_eq!(record.id, "rs1");
        assert_eq!(record.chrom, "chr1");
        assert_eq!(record.position, 101);
        assert!(record.genotype.is_call);
    }
}
