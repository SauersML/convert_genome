#![no_main]

use libfuzzer_sys::fuzz_target;
use std::io::Cursor;

fuzz_target!(|data: &[u8]| {
    // Limit records via env var for fast iteration
    std::env::set_var("CONVERT_GENOME_MAX_RECORDS", "100");

    let cursor = Cursor::new(data);
    let reader = convert_genome::dtc::Reader::new(cursor);

    // Iterate all records - should never panic
    for result in reader.take(1000) {
        match result {
            Ok(record) => {
                // Exercise Display impl
                let _ = format!("{}", record);
                // Exercise is_missing
                let _ = record.is_missing();
            }
            Err(_) => {
                // Parse errors are expected for random input
            }
        }
    }
});
