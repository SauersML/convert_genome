use noodles::vcf::variant::record::samples::keys::key as format_key;
use noodles::vcf::variant::record_buf::{Samples, samples::sample::Value};
use std::collections::HashMap;
use tracing;

/// Remap genotype indices in all samples based on mapping.
///
/// `mapping` maps Old Index -> New Index.
/// Index 0 is REF. Indices 1+ are ALTs.
pub fn remap_sample_genotypes(samples: &Samples, mapping: &HashMap<usize, usize>) -> Samples {
    let keys = samples.keys();
    let mut new_keys_vec = Vec::new();
    // Filter keys: keep GT and safe fields (GQ, DP, MIN_DP), drop allele-dependent fields (PL, AD, GP)
    // because their values order depends on allele order and simply remapping indices isn't enough for them
    // (values need permutation).
    for key in keys.as_ref().iter() {
        if key == format_key::GENOTYPE || key == "GQ" || key == "DP" || key == "MIN_DP" {
            new_keys_vec.push(key.clone());
        } else {
            tracing::debug!("Dropping field {} during polarization/normalization", key);
        }
    }

    // Create new keys object
    let new_keys: noodles::vcf::variant::record_buf::samples::Keys =
        new_keys_vec.into_iter().collect();

    let mut new_values = Vec::new();

    for sample in samples.values() {
        let mut new_sample_vals = Vec::new();

        // Iterate through valid keys and extract values from original sample
        for key in new_keys.as_ref().iter() {
            let val_opt = sample.get(key).flatten().cloned();

            if key == format_key::GENOTYPE {
                if let Some(Value::String(gt_str)) = &val_opt {
                    new_sample_vals.push(Some(Value::String(remap_gt_string(gt_str, mapping))));
                } else {
                    new_sample_vals.push(val_opt);
                }
            } else {
                new_sample_vals.push(val_opt);
            }
        }
        new_values.push(new_sample_vals);
    }

    Samples::new(new_keys, new_values)
}

/// Remap genotype indices in a GT string (e.g., "0/1" -> "1/0")
/// Handles standard separators (/, |) and ignores '.'
pub fn remap_gt_string(gt: &str, mapping: &HashMap<usize, usize>) -> String {
    let mut result = String::with_capacity(gt.len());
    let mut current_num = String::new();

    for c in gt.chars() {
        if c.is_ascii_digit() {
            current_num.push(c);
        } else {
            // End of number
            if !current_num.is_empty() {
                if let Ok(idx) = current_num.parse::<usize>() {
                    if let Some(new_idx) = mapping.get(&idx) {
                        result.push_str(&new_idx.to_string());
                    } else {
                        // Keep original if not in mapping (shouldn't happen for valid indices)
                        result.push_str(&current_num);
                    }
                } else {
                    result.push_str(&current_num);
                }
                current_num.clear();
            }
            result.push(c);
        }
    }
    // Final number
    if !current_num.is_empty() {
        if let Ok(idx) = current_num.parse::<usize>() {
            if let Some(new_idx) = mapping.get(&idx) {
                result.push_str(&new_idx.to_string());
            } else {
                result.push_str(&current_num);
            }
        } else {
            result.push_str(&current_num);
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_remap_gt_string() {
        let mut mapping = HashMap::new();
        mapping.insert(0, 1); // 0 -> 1
        mapping.insert(1, 0); // 1 -> 0
        mapping.insert(2, 2); // 2 -> 2

        assert_eq!(remap_gt_string("0/1", &mapping), "1/0");
        assert_eq!(remap_gt_string("1|2", &mapping), "0|2");
        assert_eq!(remap_gt_string("./0", &mapping), "./1");
        assert_eq!(remap_gt_string("0", &mapping), "1");
    }
}
