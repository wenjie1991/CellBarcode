use extendr_api::prelude::*;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use flate2::read::MultiGzDecoder;


// Function to read and count sequences in Fastq.gz file
#[extendr]
pub fn read_fastq_gz(in_file_path: String) -> Robj {
    // Open the gzipped file
    let file = File::open(&in_file_path).expect("Unable to open input file");
    let gz = MultiGzDecoder::new(file);
    let reader = BufReader::new(gz);


    let mut seq_map: BTreeMap<String, i32> = BTreeMap::new();

    let mut lines = reader.lines();

    while let (Some(_title), Some(seq), Some(_sep), Some(_qual)) = (
        lines.next(),
        lines.next(),
        lines.next(),
        lines.next(),
    ) {
        let seq = seq.unwrap();
        let counter = seq_map.entry(seq).or_insert(0);
        *counter += 1;
    }

    // Create frequency and sequence vectors
    let mut freq = Vec::new();
    let mut seq = Vec::new();

    for (key, value) in seq_map.iter() {
        seq.push(key.clone());
        freq.push(*value);
    }

    // Return DataFrame
    data_frame!(
        freq = freq,
        seq = seq
    )
}

// Function to read and count sequences in plain Fastq file
#[extendr]
pub fn read_fastq(in_file_path: String) -> Robj {
    // Open the plain text Fastq file
    let file = File::open(&in_file_path).expect("Unable to open input file");
    let reader = BufReader::new(file);

    let mut seq_map: BTreeMap<String, i32> = BTreeMap::new();

    let mut lines = reader.lines();
    while let (Some(_title), Some(seq), Some(_sep), Some(_qual)) = (
        lines.next(),
        lines.next(),
        lines.next(),
        lines.next(),
    ) {
        let seq = seq.unwrap();
        let counter = seq_map.entry(seq).or_insert(0);
        *counter += 1;
    }

    // Create frequency and sequence vectors
    let mut freq = Vec::new();
    let mut seq = Vec::new();

    for (key, value) in seq_map.iter() {
        seq.push(key.clone());
        freq.push(*value as i32);
    }

    // Return DataFrame
    data_frame!(
        freq = freq,
        seq = seq
    )
}

// Function to read and count sequences in paired Fastq.gz files
#[extendr]
pub fn read_fastq_gz2(in_fq1: String, in_fq2: String) -> Robj {
    // Open the gzipped files
    let file1 = File::open(&in_fq1).expect("Unable to open first input file");
    let gz1 = MultiGzDecoder::new(file1);
    let reader1 = BufReader::new(gz1);

    let file2 = File::open(&in_fq2).expect("Unable to open second input file");
    let gz2 = MultiGzDecoder::new(file2);
    let reader2 = BufReader::new(gz2);

    let mut seq_map: BTreeMap<String, i32> = BTreeMap::new();

    let mut lines1 = reader1.lines();
    let mut lines2 = reader2.lines();
    while let (
        Some(_title1), Some(seq1), Some(_), Some(_qual1),
        Some(_title2), Some(seq2), Some(_), Some(_qual2)
    ) = (
        lines1.next(),
        lines1.next(),
        lines1.next(),
        lines1.next(),
        lines2.next(),
        lines2.next(),
        lines2.next(),
        lines2.next(),
    ) {
        let seq1 = seq1.unwrap();
        let seq2 = seq2.unwrap();

        // Concatenate paired-end sequences
        let joined_seq = format!("{}{}", seq1, seq2);

        let counter = seq_map.entry(joined_seq).or_insert(0);
        *counter += 1;
    }

    // Create frequency and sequence vectors
    let mut freq = Vec::new();
    let mut seq = Vec::new();

    for (key, value) in seq_map.iter() {
        seq.push(key.clone());
        freq.push(*value);
    }

    // Return DataFrame
    data_frame!(
        freq = freq,
        seq = seq
    )
}

extendr_module! {
    mod lib_read_seq;
    fn read_fastq;
    fn read_fastq_gz;
    fn read_fastq_gz2;
}
