use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;
use std::fs::File;
use flate2::read::{MultiGzDecoder, DeflateDecoder};

/// Opens a file and transparently peels off GZIP, BGZF, and ZIP layers
/// to expose the underlying raw data stream.
///
/// Supports nested compression (e.g. .vcf.gz.zip).
/// For ZIP files, it reads the first entry.
pub fn open_input(path: &Path) -> anyhow::Result<Box<dyn BufRead + Send>> {
    let file = File::open(path)?;
    let mut reader: Box<dyn BufRead + Send> = Box::new(BufReader::new(file));

    // Limit recursion depth to avoid infinite loops on malformed inputs
    let mut depth = 0;
    const MAX_DEPTH: usize = 10;

    loop {
        if depth >= MAX_DEPTH {
            break;
        }

        let is_gzip;
        let is_zip;
        {
            let buf = reader.fill_buf()?;
            if buf.is_empty() {
                break;
            }
            // GZIP magic: 1f 8b
            is_gzip = buf.len() >= 2 && buf[0] == 0x1f && buf[1] == 0x8b;
            // ZIP magic: PK\x03\x04 (Little Endian: 50 4b 03 04)
            is_zip = buf.len() >= 4 && buf[0] == 0x50 && buf[1] == 0x4b && buf[2] == 0x03 && buf[3] == 0x04;
        }

        if is_gzip {
            tracing::debug!("Detected GZIP/BGZF layer");
            // Use MultiGzDecoder to support BGZF and concatenated GZIP members
            reader = Box::new(BufReader::new(MultiGzDecoder::new(reader)));
            depth += 1;
        } else if is_zip {
            tracing::debug!("Detected ZIP layer");
            reader = Box::new(BufReader::new(read_zip_first_entry(reader)?));
            depth += 1;
        } else {
            // No known compression found
            break;
        }
    }
    Ok(reader)
}

/// Reads the first entry of a ZIP stream.
fn read_zip_first_entry(mut reader: Box<dyn BufRead + Send>) -> anyhow::Result<Box<dyn Read + Send>> {
    // Read 30 bytes Local File Header
    let mut header = [0u8; 30];
    reader.read_exact(&mut header)?;

    // Parse fields
    // Signature (0-3) should be 50 4b 03 04, detected by caller
    let flags = u16::from_le_bytes([header[6], header[7]]);
    let compression = u16::from_le_bytes([header[8], header[9]]);
    let compressed_size = u32::from_le_bytes([header[18], header[19], header[20], header[21]]) as u64;
    // let uncompressed_size = u32::from_le_bytes([header[22], header[23], header[24], header[25]]);
    let name_len = u16::from_le_bytes([header[26], header[27]]) as usize;
    let extra_len = u16::from_le_bytes([header[28], header[29]]) as usize;

    // Skip filename and extra field
    // We must read and discard since we can't seek
    if name_len > 0 {
        io::copy(&mut reader.by_ref().take(name_len as u64), &mut io::sink())?;
    }
    if extra_len > 0 {
        io::copy(&mut reader.by_ref().take(extra_len as u64), &mut io::sink())?;
    }

    match compression {
        8 => {
            // Deflate
            // DeflateDecoder reads until end of deflate block.
            // It consumes the compressed data from `reader`.
            Ok(Box::new(DeflateDecoder::new(reader)))
        }
        0 => {
            // Store (No compression)
            // If bit 3 (0x0008) of flags is set, size is in Data Descriptor (after data).
            // For Store method, this is problematic for streaming if we don't know when to stop.
            if (flags & 0x0008) != 0 {
                tracing::warn!("ZIP Stored entry has Data Descriptor (bit 3 set). Cannot determine size safely in stream. Reading until EOF.");
                Ok(Box::new(reader))
            } else {
                // Size is known in header
                Ok(Box::new(reader.take(compressed_size)))
            }
        }
        _ => {
            anyhow::bail!("Unsupported ZIP compression method: {}", compression);
        }
    }
}
