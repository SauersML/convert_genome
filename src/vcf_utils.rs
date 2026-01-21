use noodles::vcf::variant::record::samples::keys::key as format_key;
use noodles::vcf::variant::record_buf::{Samples, samples::sample::Value};
use std::collections::HashMap;
use tracing;

/// Remap genotype indices in all samples based on mapping.
///
/// `mapping` maps Old Index -> New Index.
/// Index 0 is REF. Indices 1+ are ALTs.
/// Remap genotype indices in all samples based on mapping.
///
/// `mapping` maps Old Index -> New Index.
/// Index 0 is REF. Indices 1+ are ALTs.
pub fn remap_sample_genotypes(samples: &Samples, mapping: &HashMap<usize, usize>) -> Samples {
    let keys = samples.keys();
    
    // First, verify we have keys and map the fields we care about to their NUMERIC INDICES in the input
    let mut gt_index = None;
    let mut other_indices = Vec::new();
    
    // Create the new keys list
    let mut new_keys_vec = Vec::new();
    
    for (i, key) in keys.as_ref().iter().enumerate() {
        if key == format_key::GENOTYPE {
            gt_index = Some(i);
            new_keys_vec.push(key.clone());
        } else if key == "GQ" || key == "DP" || key == "MIN_DP" {
            other_indices.push((i, key.clone()));
            new_keys_vec.push(key.clone());
        } else {
            // Drop other fields (PL, AD, GP) as they become invalid when alleles change
            // and we cannot easily reorder/recalculate them without more complex logic.
            tracing::debug!("Dropping field {} during polarization/normalization", key);
        }
    }

    // Create new keys object
    let new_keys: noodles::vcf::variant::record_buf::samples::Keys =
        new_keys_vec.into_iter().collect();

    let mut new_values = Vec::new();

    for sample in samples.values() {
        let mut new_sample_vals = Vec::new();

        // 1. Process GT if present (always first in standard VCFs, but let's follow our key order)
        if let Some(idx) = gt_index {
            // Look up value by NUMERIC INDEX (idx) using direct slice access
            let val_opt = sample.values().get(idx).cloned().flatten();
            
            if let Some(Value::String(gt_str)) = &val_opt {
                let new_gt = remap_gt_string(gt_str, mapping);
                if !mapping.is_empty() && gt_str != &new_gt {
                    tracing::debug!("Remapped GT: {} -> {} (mapping: {:?})", gt_str, new_gt, mapping);
                }
                new_sample_vals.push(Some(Value::String(new_gt)));
            } else {
                new_sample_vals.push(val_opt);
            }
        }
        
        // 2. Process other preserved fields
        for (idx, _) in &other_indices {
             let val_opt = sample.values().get(*idx).cloned().flatten();
             new_sample_vals.push(val_opt);
        }
        
        new_values.push(new_sample_vals);
    }

    Samples::new(new_keys, new_values)
}

/// Remap genotype indices in a GT string (e.g., "0/1" -> "1/0")
/// Handles standard separators (/, |) and ignores '.'
pub fn remap_gt_string(gt: &str, mapping: &HashMap<usize, usize>) -> String {
    // Fast path: if all allele indices are single ASCII digits (0-9), remap byte-by-byte.
    // This covers the overwhelmingly common cases like 0/0, 0/1, 1|0, ./., ./0, etc.
    // If we detect any multi-digit index, fall back to the general parser below.
    if gt
        .as_bytes()
        .iter()
        .all(|&b| matches!(b, b'0'..=b'9' | b'.' | b'/' | b'|'))
    {
        // Detect multi-digit indices by checking for adjacent digits.
        let bytes = gt.as_bytes();
        let mut has_multi_digit = false;
        for w in bytes.windows(2) {
            if w[0].is_ascii_digit() && w[1].is_ascii_digit() {
                has_multi_digit = true;
                break;
            }
        }

        if !has_multi_digit {
            let mut out: Vec<u8> = Vec::with_capacity(bytes.len());
            for &b in bytes {
                if b.is_ascii_digit() {
                    let idx = (b - b'0') as usize;
                    if let Some(&new_idx) = mapping.get(&idx)
                        && new_idx < 10
                    {
                        out.push((new_idx as u8) + b'0');
                    } else {
                        // Keep original if not in mapping or not representable as single digit.
                        out.push(b);
                    }
                } else {
                    out.push(b);
                }
            }
            return String::from_utf8(out).unwrap_or_else(|_| gt.to_string());
        }
    }

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
