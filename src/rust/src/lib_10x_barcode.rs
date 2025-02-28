use extendr_api::prelude::*;
use regex::Regex;
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader};

// Define Barcode struct
#[derive(Debug, Clone)]
struct Barcode {
    cell_barcode: String,
    umi: String,
    barcode: String,
    count: usize,
}

/// Parse 10X SAM file
/// 
/// @param in_file_path A string, define the un-mapped sequences
/// @param regex_str A string, define the regular expression to match the barcode sequence.
/// @param cell_barcode_tag A string, define the tag of 10X cell barcode field in sam file.
/// @param umi_tag A string, define the tag of UMI field in the sam file.
/// @return A list containing two data frames:
///   - `barcode_df`: with columns `cell_barcode`, `umi`, `barcode_seq`, and `count`.
///   - `raw_reads_df`: with columns `cell_barcode` and `count`.
#[extendr]
pub fn parse_10x_sam(
    in_file_path: String,
    regex_str: String,
    cell_barcode_tag: String,
    umi_tag: String,
) -> Robj {
    // Compile the regular expression
    let str_expr = Regex::new(&regex_str).expect("Invalid regular expression.");

    // Open the input file
    let infile = File::open(&in_file_path).expect("Unable to open input file.");
    let reader = BufReader::new(infile);

    // Data containers
    let mut seq_map: BTreeMap<String, Barcode> = BTreeMap::new();
    let mut raw_count: BTreeMap<String, usize> = BTreeMap::new();

    // Parse each line
    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with('@') {
            continue;
        }
        
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 12 {
            continue;
        }

        let seq = parts[9];
        let mut cell_barcode = String::new();
        let mut umi_seq = String::new();
        
        // Extract cell barcode and UMI
        for part in &parts[11..] {
            if part.starts_with(&cell_barcode_tag) {
                cell_barcode = part[5..].to_string();
            } else if part.starts_with(&umi_tag) {
                umi_seq = part[5..].to_string();
            }
        }

        // Count raw reads
        let count = raw_count.entry(cell_barcode.clone()).or_insert(0);
        *count += 1;

        // Match barcode using regex
        if let Some(caps) = str_expr.captures(seq) {
            if let Some(matched) = caps.get(1) {
                let barcode = matched.as_str().to_string();
                let key = format!("{}{}{}", cell_barcode, umi_seq, barcode);

                // Store in seq_map
                seq_map.entry(key.clone())
                    .and_modify(|e| e.count += 1)
                    .or_insert(Barcode {
                        cell_barcode: cell_barcode.clone(),
                        umi: umi_seq.clone(),
                        barcode: barcode.clone(),
                        count: 1,
                    });
            }
        }
    }

    // Create DataFrames for output

    // Barcode DataFrame
    let mut cell_barcodes = Vec::new();
    let mut umis = Vec::new();
    let mut barcodes = Vec::new();
    let mut counts = Vec::new();
    let mut unique_cell_barcodes = BTreeSet::new();

    for barcode in seq_map.values() {
        cell_barcodes.push(barcode.cell_barcode.clone());
        umis.push(barcode.umi.clone());
        barcodes.push(barcode.barcode.clone());
        counts.push(barcode.count as i32);
        unique_cell_barcodes.insert(barcode.cell_barcode.clone());
    }

    let barcode_df = data_frame!( 
        cell_barcode = cell_barcodes, 
        umi = umis, 
        barcode_seq = barcodes, 
        count = counts 
    );

    // Raw Reads DataFrame
    let mut unique_barcodes: Vec<String> = Vec::new();
    let mut raw_counts: Vec<i32> = Vec::new();

    for cell_barcode in unique_cell_barcodes {
        unique_barcodes.push(cell_barcode.clone());
        raw_counts.push(*raw_count.get(&cell_barcode).unwrap_or(&0) as i32);
    }

    let raw_reads_df = data_frame!(
        cell_barcode = unique_barcodes,
        count = raw_counts
    );

    // Return list of data frames
    list!(barcode_df = barcode_df, raw_reads_df = raw_reads_df).into()
}

extendr_module! {
    mod lib_10x_barcode;
    fn parse_10x_sam;
}

