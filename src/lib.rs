#![doc = include_str!("../README.md")]

pub mod cli;
pub mod conversion;
pub mod dtc;
pub mod reference;

pub use conversion::{ConversionConfig, ConversionSummary, OutputFormat, convert_dtc_file};
