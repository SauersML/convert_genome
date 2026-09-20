//! Placing one standardized record on the reference panel's alleles, and writing
//! the panel sites a complete call set left out.
//!
//! Two facts about the input steer this, and both travel in the input's own VCF
//! header rather than on the command line, because they describe the data:
//!
//! * `##pgsStrandFlipPrior=<p>` — the probability that any one site of this
//!   input is reported on the opposite strand. It is 0 for a sequencing VCF
//!   (strand fixed by alignment) and is measured from the input's own
//!   non-palindromic heterozygotes otherwise. It decides whether an
//!   unambiguous flip is complemented and how A/T and C/G homozygotes are
//!   resolved (see [`crate::harmonize::resolve_snv`]). Absent: 0.
//! * `##pgsAbsentGenotypes=HomRef` — the input is a complete call set (a
//!   high-coverage whole-genome VCF): a site it does not mention is homozygous
//!   reference, except inside a no-call block (a gVCF-style record with END and
//!   a missing genotype) or under a called deletion. Every panel record the
//!   input did not write is then emitted, so an imputation run downstream has
//!   the genome's own genotype at every panel site and imputes only true
//!   no-calls. Absent: nothing is filled, which is right for an array.

use anyhow::Result;
use noodles::core::Position;
use noodles::vcf::{
    self,
    header::record::value::{Collection, map::info::Number},
    variant::record::samples::series::value::genotype::Phasing,
    variant::record_buf::{AlternateBases, RecordBuf, Samples, samples::sample::Value},
};

use crate::harmonize::{Resolution, complement_seq, resolve_snv};
use crate::panel::{ChromClaims, PaddedPanel, PanelSite};

/// Header key: `HomRef` means sites absent from the input are homozygous reference.
pub const ABSENT_GENOTYPES_KEY: &str = "pgsAbsentGenotypes";
/// Header key: the per-site strand-flip prior of the input.
pub const STRAND_FLIP_PRIOR_KEY: &str = "pgsStrandFlipPrior";

/// What the input header says about the input.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct InputQcMeta {
    pub absent_hom_ref: bool,
    pub strand_prior: f64,
}

impl InputQcMeta {
    pub fn from_header(header: &vcf::Header) -> Self {
        let value = |key: &str| -> Option<String> {
            match header.other_records().get(key)? {
                Collection::Unstructured(values) => values.last().cloned(),
                Collection::Structured(_) => None,
            }
        };
        let absent_hom_ref =
            value(ABSENT_GENOTYPES_KEY).is_some_and(|v| v.trim().eq_ignore_ascii_case("HomRef"));
        let strand_prior = value(STRAND_FLIP_PRIOR_KEY)
            .and_then(|v| v.trim().parse::<f64>().ok())
            .filter(|p| p.is_finite() && *p >= 0.0)
            .map(|p| p.min(1.0))
            .unwrap_or(0.0);
        Self {
            absent_hom_ref,
            strand_prior,
        }
    }
}

/// The first sample's genotype: allele indices (`None` = missing) and the
/// separators between them.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Gt {
    alleles: Vec<Option<usize>>,
    separators: Vec<char>,
}

impl Gt {
    fn missing_like(&self) -> Self {
        Self {
            alleles: vec![None; self.alleles.len()],
            separators: self.separators.clone(),
        }
    }

    fn render(&self) -> String {
        let mut out = String::new();
        for (i, allele) in self.alleles.iter().enumerate() {
            if i > 0 {
                out.push(self.separators.get(i - 1).copied().unwrap_or('/'));
            }
            match allele {
                Some(n) => out.push_str(&n.to_string()),
                None => out.push('.'),
            }
        }
        out
    }
}

fn parse_gt_str(s: &str) -> Option<Gt> {
    let mut alleles = Vec::new();
    let mut separators = Vec::new();
    let mut token = String::new();
    for c in s.chars() {
        if c == '/' || c == '|' {
            alleles.push(parse_allele(&token)?);
            separators.push(c);
            token.clear();
        } else {
            token.push(c);
        }
    }
    alleles.push(parse_allele(&token)?);
    Some(Gt {
        alleles,
        separators,
    })
}

fn parse_allele(token: &str) -> Option<Option<usize>> {
    if token == "." {
        return Some(None);
    }
    token.parse::<usize>().ok().map(Some)
}

fn first_sample_gt(record: &RecordBuf) -> Option<Gt> {
    let samples = record.samples();
    let gt_index = samples.keys().as_ref().get_index_of("GT")?;
    let sample = samples.values().next()?;
    match sample.values().get(gt_index)?.as_ref()? {
        Value::String(s) => parse_gt_str(s),
        Value::Genotype(genotype) => {
            let mut alleles = Vec::new();
            let mut separators = Vec::new();
            for (i, allele) in genotype.as_ref().iter().enumerate() {
                if i > 0 {
                    separators.push(match allele.phasing() {
                        Phasing::Phased => '|',
                        Phasing::Unphased => '/',
                    });
                }
                alleles.push(allele.position());
            }
            (!alleles.is_empty()).then_some(Gt {
                alleles,
                separators,
            })
        }
        _ => None,
    }
}

/// The first sample's fields with GT replaced; GQ/DP/MIN_DP survive, every
/// allele-indexed field (AD, PL, GP) is dropped because the alleles changed.
fn samples_with_gt(samples: &Samples, gt: String) -> Samples {
    let keys = samples.keys().as_ref();
    let kept: Vec<(usize, &String)> = keys
        .iter()
        .enumerate()
        .filter(|(_, k)| matches!(k.as_str(), "GQ" | "DP" | "MIN_DP"))
        .collect();
    let first = samples.values().next();
    let mut values = vec![Some(Value::String(gt))];
    for (i, _) in &kept {
        values.push(
            first
                .as_ref()
                .and_then(|s| s.values().get(*i).cloned().flatten()),
        );
    }
    let new_keys = std::iter::once("GT".to_string())
        .chain(kept.iter().map(|(_, k)| (*k).clone()))
        .collect();
    Samples::new(new_keys, vec![values])
}

fn is_symbolic_block_allele(alt: &str) -> bool {
    alt == "<*>" || alt.eq_ignore_ascii_case("<NON_REF>") || alt == "*" || alt == "."
}

fn info_end(record: &RecordBuf) -> Option<u64> {
    use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
    match record.info().get("END")? {
        Some(InfoValue::Integer(n)) => u64::try_from(*n).ok(),
        _ => None,
    }
}

/// Handle a gVCF-style block record (END in INFO, no real ALT) in fill mode.
///
/// A block with a missing genotype becomes a no-call span: panel sites inside
/// it are written as missing, never as homozygous reference. A block with a
/// called genotype is reference-confident and needs nothing (absent sites are
/// already homozygous reference). Either way the block itself is not emitted.
/// Returns true when `record` was such a block.
pub fn take_block_record(record: &RecordBuf, panel: &mut PaddedPanel, output_name: &str) -> bool {
    let alts = record.alternate_bases().as_ref();
    if !alts.iter().all(|a| is_symbolic_block_allele(a)) {
        return false;
    }
    let Some(end) = info_end(record) else {
        return false;
    };
    let Some(start) = record.variant_start().map(|p| usize::from(p) as u64) else {
        return false;
    };
    let no_call = first_sample_gt(record).is_none_or(|gt| gt.alleles.iter().all(Option::is_none));
    if no_call {
        let chrom = record.reference_sequence_name().to_string();
        if let Some(claims) = panel.claims_for(&chrom, output_name) {
            claims.no_call_spans.push((start, end.max(start)));
        }
    }
    true
}

/// Bases a called non-SNV variant changes: the whole REF span, or all but the
/// shared anchor base when every ALT starts with it (a left-anchored indel).
fn altered_span(pos: u64, reference: &str, alts: &[String]) -> Option<(u64, u64)> {
    let len = reference.len() as u64;
    let anchored = !alts.is_empty()
        && alts
            .iter()
            .all(|a| !a.is_empty() && a[..1].eq_ignore_ascii_case(&reference[..1]));
    let start = if anchored { pos + 1 } else { pos };
    let end = pos + len.saturating_sub(1);
    (start <= end).then_some((start, end))
}

fn claims<'a>(
    panel: &'a mut PaddedPanel,
    meta: &InputQcMeta,
    chrom: &str,
) -> Option<&'a mut ChromClaims> {
    if !meta.absent_hom_ref {
        return None;
    }
    panel.claims_for(chrom, chrom)
}

fn strip_allele_indexed_info(
    record: &RecordBuf,
    header: &vcf::Header,
) -> vcf::variant::record_buf::Info {
    let mut info = record.info().clone();
    let infos = header.infos();
    info.as_mut().retain(|key, _| match infos.get(key) {
        Some(definition) => !matches!(
            definition.number(),
            Number::AlternateBases | Number::ReferenceAlternateBases | Number::Samples
        ),
        None => true,
    });
    info
}

fn rebuild(
    record: &RecordBuf,
    header: &vcf::Header,
    reference: &str,
    alts: &[String],
    gt: String,
) -> RecordBuf {
    let mut builder = RecordBuf::builder()
        .set_reference_sequence_name(record.reference_sequence_name())
        .set_ids(record.ids().clone())
        .set_filters(record.filters().clone())
        .set_reference_bases(reference.to_string())
        .set_alternate_bases(AlternateBases::from(alts.to_vec()))
        .set_info(strip_allele_indexed_info(record, header))
        .set_samples(samples_with_gt(record.samples(), gt));
    if let Some(start) = record.variant_start() {
        builder = builder.set_variant_start(start);
    }
    if let Some(qual) = record.quality_score() {
        builder = builder.set_quality_score(qual);
    }
    builder.build()
}

fn with_missing_gt(record: RecordBuf, gt: &Gt) -> RecordBuf {
    let mut record = record;
    let samples = samples_with_gt(record.samples(), gt.missing_like().render());
    *record.samples_mut() = samples;
    record
}

/// Place one standardized record on the panel's alleles.
///
/// Returns the record to emit. SNVs at panel positions are resolved per site
/// (strand, palindromes, split records); a call that cannot be placed is
/// emitted with a missing genotype, never on a guessed allele. Indels and
/// sites the panel lacks pass through unchanged. `track_padding` records novel
/// alleles and sites for the padded-panel output; it is off in the pipeline,
/// where registering every private variant cost O(n^2).
pub fn apply_panel(
    record: RecordBuf,
    panel: &mut PaddedPanel,
    header: &vcf::Header,
    meta: &InputQcMeta,
    track_padding: bool,
) -> RecordBuf {
    let chrom = record.reference_sequence_name().to_string();
    let Some(pos) = record.variant_start().map(|p| usize::from(p) as u64) else {
        return record;
    };
    let record_ref = record.reference_bases().to_ascii_uppercase();
    let record_alts: Vec<String> = record
        .alternate_bases()
        .as_ref()
        .iter()
        .map(|a| a.to_ascii_uppercase())
        .collect();
    let gt = first_sample_gt(&record);
    let allele_text = |index: usize| -> Option<String> {
        if index == 0 {
            Some(record_ref.clone())
        } else {
            record_alts.get(index - 1).cloned()
        }
    };
    let called: Option<Vec<String>> = gt.as_ref().and_then(|gt| {
        gt.alleles
            .iter()
            .map(|a| a.and_then(allele_text))
            .collect::<Option<Vec<String>>>()
    });
    let partly_called: Vec<String> = gt
        .as_ref()
        .map(|gt| {
            gt.alleles
                .iter()
                .filter_map(|a| a.and_then(allele_text))
                .collect()
        })
        .unwrap_or_default();
    let is_snv = record_ref.len() == 1 && record_alts.iter().all(|a| a.len() == 1);
    let panel_records: Vec<PanelSite> = panel.original().get_all(&chrom, pos).to_vec();

    if panel_records.is_empty() {
        if track_padding {
            for allele in std::iter::once(&record_ref).chain(record_alts.iter()) {
                panel.get_or_add_allele_index(&chrom, pos, allele, &record_ref);
            }
        }
        panel.qc.absent_from_panel += 1;
        claim_unplaced(
            panel,
            meta,
            &chrom,
            pos,
            &record_ref,
            &record_alts,
            is_snv,
            called,
        );
        return record;
    }

    if !is_snv {
        let exact = panel_records.iter().find(|r| {
            r.ref_allele.eq_ignore_ascii_case(&record_ref)
                && r.alt_alleles.len() == record_alts.len()
                && r.alt_alleles
                    .iter()
                    .zip(&record_alts)
                    .all(|(a, b)| a.eq_ignore_ascii_case(b))
        });
        if exact.is_some() {
            panel.qc.indel_exact += 1;
        } else {
            panel.qc.indel_unmatched += 1;
        }
        let exact_key = exact.map(|r| emitted_key(pos, &r.ref_allele, &r.alt_alleles));
        if let Some(c) = claims(panel, meta, &chrom)
            && let Some(key) = exact_key
        {
            c.emitted.insert(key);
        }
        claim_unplaced(
            panel,
            meta,
            &chrom,
            pos,
            &record_ref,
            &record_alts,
            false,
            called,
        );
        return record;
    }

    let candidates: Vec<&PanelSite> = panel_records
        .iter()
        .filter(|r| r.is_snv() && r.ref_allele.eq_ignore_ascii_case(&record_ref))
        .collect();
    let Some(gt) = gt else {
        return record;
    };

    if candidates.is_empty() {
        if panel_records.iter().any(PanelSite::is_snv) {
            // The panel's REF at this position is not the record's (which is
            // the reference FASTA's after standardization). No safe mapping.
            panel.qc.snv_reference_mismatch += 1;
            if let Some(c) = claims(panel, meta, &chrom) {
                c.snv_calls.insert(pos, None);
            }
            return with_missing_gt(record, &gt);
        }
        panel.qc.absent_from_panel += 1;
        claim_unplaced(
            panel,
            meta,
            &chrom,
            pos,
            &record_ref,
            &record_alts,
            true,
            called,
        );
        return record;
    }

    if partly_called.is_empty() {
        panel.qc.no_call += 1;
        return missing_on_panel_alleles(
            record,
            panel,
            header,
            meta,
            &chrom,
            pos,
            &gt,
            &candidates,
        );
    }

    let mut union: Vec<String> = Vec::new();
    for site in &candidates {
        for alt in &site.alt_alleles {
            let alt = alt.to_ascii_uppercase();
            if !union.contains(&alt) {
                union.push(alt);
            }
        }
    }
    let alt_frequency = if union.len() == 1 {
        candidates.iter().find_map(|site| {
            let idx = site
                .alt_alleles
                .iter()
                .position(|a| a.eq_ignore_ascii_case(&union[0]))?;
            site.alt_frequency(idx)
        })
    } else {
        None
    };
    let (resolution, class) = resolve_snv(
        &partly_called,
        &record_ref,
        &union,
        alt_frequency,
        meta.strand_prior,
    );
    panel.qc.count(class);

    let resolved: Vec<Option<String>> = match resolution {
        Resolution::Missing => {
            return missing_on_panel_alleles(
                record,
                panel,
                header,
                meta,
                &chrom,
                pos,
                &gt,
                &candidates,
            );
        }
        Resolution::Keep => gt.alleles.iter().map(|a| a.and_then(allele_text)).collect(),
        Resolution::Complement => gt
            .alleles
            .iter()
            .map(|a| a.and_then(allele_text).map(|s| complement_seq(&s)))
            .collect(),
    };

    let mut non_ref: Vec<&String> = Vec::new();
    for allele in resolved.iter().flatten() {
        if !allele.eq_ignore_ascii_case(&record_ref) && !non_ref.contains(&allele) {
            non_ref.push(allele);
        }
    }
    let chosen = if non_ref.is_empty() {
        // A homozygous-reference call: place it on the most common panel ALT,
        // the one an array probe at this position was most likely designed for.
        candidates.iter().copied().max_by(|a, b| {
            let fa = a.alt_frequency(0).unwrap_or(-1.0);
            let fb = b.alt_frequency(0).unwrap_or(-1.0);
            fa.partial_cmp(&fb)
                .unwrap_or(std::cmp::Ordering::Equal)
                // max_by keeps the last maximum; prefer the first record on ties.
                .then(std::cmp::Ordering::Greater)
        })
    } else {
        candidates.iter().copied().find(|site| {
            non_ref
                .iter()
                .all(|n| site.alt_alleles.iter().any(|a| a.eq_ignore_ascii_case(n)))
        })
    };
    let Some(site) = chosen else {
        // Two non-reference alleles that live on different split records: no
        // single panel record can carry the genotype.
        panel.qc.split_record_het += 1;
        return missing_on_panel_alleles(
            record,
            panel,
            header,
            meta,
            &chrom,
            pos,
            &gt,
            &candidates,
        );
    };

    let site_ref = site.ref_allele.to_ascii_uppercase();
    let site_alts: Vec<String> = site
        .alt_alleles
        .iter()
        .map(|a| a.to_ascii_uppercase())
        .collect();
    let new_gt = Gt {
        alleles: resolved
            .iter()
            .map(|a| {
                a.as_ref().map(|allele| {
                    if *allele == site_ref {
                        0
                    } else {
                        site_alts
                            .iter()
                            .position(|x| x == allele)
                            .map_or(0, |i| i + 1)
                    }
                })
            })
            .collect(),
        separators: gt.separators.clone(),
    };
    let emitted = emitted_key(pos, &site_ref, &site_alts);
    let called_after: Option<Vec<String>> = resolved.iter().cloned().collect();
    if let Some(c) = claims(panel, meta, &chrom) {
        c.emitted.insert(emitted);
        c.snv_calls.insert(pos, called_after);
    }
    rebuild(&record, header, &site_ref, &site_alts, new_gt.render())
}

/// Emit a genotype that could not be placed as missing, on the first panel
/// record's alleles, so an imputation run treats the position as a known panel
/// marker to impute rather than an unknown one to discard.
#[allow(clippy::too_many_arguments)]
fn missing_on_panel_alleles(
    record: RecordBuf,
    panel: &mut PaddedPanel,
    header: &vcf::Header,
    meta: &InputQcMeta,
    chrom: &str,
    pos: u64,
    gt: &Gt,
    candidates: &[&PanelSite],
) -> RecordBuf {
    let Some(site) = candidates.first() else {
        if let Some(c) = claims(panel, meta, chrom) {
            c.snv_calls.insert(pos, None);
        }
        return with_missing_gt(record, gt);
    };
    let site_ref = site.ref_allele.to_ascii_uppercase();
    let site_alts: Vec<String> = site
        .alt_alleles
        .iter()
        .map(|a| a.to_ascii_uppercase())
        .collect();
    let key = emitted_key(pos, &site_ref, &site_alts);
    if let Some(c) = claims(panel, meta, chrom) {
        c.emitted.insert(key);
        c.snv_calls.insert(pos, None);
    }
    rebuild(
        &record,
        header,
        &site_ref,
        &site_alts,
        gt.missing_like().render(),
    )
}

fn emitted_key(pos: u64, reference: &str, alts: &[String]) -> (u64, String, String) {
    (
        pos,
        reference.to_ascii_uppercase(),
        alts.iter()
            .map(|a| a.to_ascii_uppercase())
            .collect::<Vec<_>>()
            .join(","),
    )
}

/// Record what a record the panel cannot place says about its position.
#[allow(clippy::too_many_arguments)]
fn claim_unplaced(
    panel: &mut PaddedPanel,
    meta: &InputQcMeta,
    chrom: &str,
    pos: u64,
    reference: &str,
    alts: &[String],
    is_snv: bool,
    called: Option<Vec<String>>,
) {
    let Some(c) = claims(panel, meta, chrom) else {
        return;
    };
    if is_snv {
        c.snv_calls.insert(pos, called);
        return;
    }
    let carries_variant = match &called {
        Some(alleles) => alleles.iter().any(|a| !a.eq_ignore_ascii_case(reference)),
        None => true,
    };
    if carries_variant && let Some(span) = altered_span(pos, reference, alts) {
        c.altered_spans.push(span);
    }
}

fn merge_spans(mut spans: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    spans.sort_unstable();
    let mut merged: Vec<(u64, u64)> = Vec::with_capacity(spans.len());
    for (start, end) in spans {
        match merged.last_mut() {
            Some(last) if start <= last.1.saturating_add(1) => last.1 = last.1.max(end),
            _ => merged.push((start, end)),
        }
    }
    merged
}

fn overlaps(spans: &[(u64, u64)], start: u64, end: u64) -> bool {
    let i = spans.partition_point(|&(_, e)| e < start);
    spans.get(i).is_some_and(|&(s, _)| s <= end)
}

/// Emit every panel record the input did not write, for each chromosome the
/// input touched, as homozygous reference (or missing, inside a no-call block
/// or under a called deletion, or at a position whose own call was missing).
pub fn fill_absent_panel_records(
    panel: &mut PaddedPanel,
    mut push: impl FnMut(RecordBuf) -> Result<()>,
) -> Result<()> {
    let all_claims = std::mem::take(&mut panel.claims);
    let (index, qc) = panel.original_and_qc();
    let gt_keys: noodles::vcf::variant::record_buf::samples::Keys =
        std::iter::once("GT".to_string()).collect();
    for (panel_chrom, claims) in all_claims {
        let no_call = merge_spans(claims.no_call_spans);
        let altered = merge_spans(claims.altered_spans);
        for &pos in index.positions(&panel_chrom) {
            for site in index.get_all(&panel_chrom, pos) {
                let key = emitted_key(pos, &site.ref_allele, &site.alt_alleles);
                if claims.emitted.contains(&key) {
                    continue;
                }
                let end = pos + (site.ref_allele.len() as u64).saturating_sub(1);
                let gt = if overlaps(&no_call, pos, end) {
                    qc.fill_missing_no_call += 1;
                    "./.".to_string()
                } else if overlaps(&altered, pos, end) {
                    qc.fill_missing_overlap += 1;
                    "./.".to_string()
                } else {
                    match claims.snv_calls.get(&pos) {
                        Some(None) => {
                            qc.fill_missing_no_call += 1;
                            "./.".to_string()
                        }
                        Some(Some(alleles)) if site.is_snv() => {
                            // Split representation: an allele that is not this
                            // record's ALT counts as its REF.
                            let codes: Vec<String> = alleles
                                .iter()
                                .map(|a| {
                                    site.alt_alleles
                                        .iter()
                                        .position(|x| x.eq_ignore_ascii_case(a))
                                        .map_or(0, |i| i + 1)
                                        .to_string()
                                })
                                .collect();
                            let gt = codes.join("/");
                            if codes.iter().all(|c| c == "0") {
                                qc.fill_hom_ref += 1;
                            }
                            gt
                        }
                        _ => {
                            qc.fill_hom_ref += 1;
                            "0/0".to_string()
                        }
                    }
                };
                let Ok(start) = Position::try_from(pos as usize) else {
                    continue;
                };
                let record = RecordBuf::builder()
                    .set_reference_sequence_name(claims.output_name.clone())
                    .set_variant_start(start)
                    .set_reference_bases(site.ref_allele.clone())
                    .set_alternate_bases(AlternateBases::from(site.alt_alleles.clone()))
                    .set_samples(Samples::new(
                        gt_keys.clone(),
                        vec![vec![Some(Value::String(gt))]],
                    ))
                    .build();
                push(record)?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn genotype_strings_round_trip() {
        let gt = parse_gt_str("0|1").unwrap();
        assert_eq!(gt.alleles, vec![Some(0), Some(1)]);
        assert_eq!(gt.render(), "0|1");
        assert_eq!(gt.missing_like().render(), ".|.");
        assert_eq!(parse_gt_str("1").unwrap().render(), "1");
        assert_eq!(parse_gt_str("./1").unwrap().alleles, vec![None, Some(1)]);
        assert!(parse_gt_str("A/G").is_none());
    }

    #[test]
    fn deletion_alters_everything_but_its_anchor() {
        assert_eq!(
            altered_span(100, "ACGT", &["A".to_string()]),
            Some((101, 103))
        );
        assert_eq!(altered_span(100, "A", &["AT".to_string()]), None);
        assert_eq!(
            altered_span(100, "AC", &["GT".to_string()]),
            Some((100, 101))
        );
    }

    #[test]
    fn span_lookup() {
        let spans = merge_spans(vec![(10, 20), (15, 30), (40, 40), (31, 35)]);
        assert_eq!(spans, vec![(10, 35), (40, 40)]);
        assert!(overlaps(&spans, 35, 36));
        assert!(!overlaps(&spans, 36, 39));
        assert!(overlaps(&spans, 38, 45));
        assert!(!overlaps(&spans, 1, 9));
    }
}
