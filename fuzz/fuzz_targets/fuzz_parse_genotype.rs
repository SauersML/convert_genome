#![no_main]

use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    // Convert bytes to lossy UTF-8 string for genotype parsing
    let input = String::from_utf8_lossy(data);

    // Fuzz the genotype parser - should never panic
    let alleles = convert_genome::dtc::parse_genotype(&input);

    // Invariants: result should be bounded and consistent
    assert!(alleles.len() <= input.len() + 1, "allele count explosion");

    // Each allele should be representable
    for allele in &alleles {
        let _ = format!("{:?}", allele);
    }
});
