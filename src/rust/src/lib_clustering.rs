use extendr_api::prelude::*;
use leven_distance::Levenshtein;
//use std::collections::HashMap;

fn crates_levenshtein_distance(
    source: &str,
    target: &str,
    lev: &Levenshtein,
) -> i32 {
    //let mut lev = Levenshtein::new();
    //lev.set_insert_cost(insert_cost);
    //lev.set_delete_cost(delete_cost);
    //lev.set_replace_cost(replace_cost);

    lev.calculate(source, target).distance() 
}

// Hamming Distance Function
fn hamm_dist(s1: &str, s2: &str) -> i32 {
    if s1.len() != s2.len() {
        return i32::MAX;
    }
    s1.chars().zip(s2.chars()).filter(|(c1, c2)| c1 != c2).count() as i32
}


/// Sequence clustering using Levenshtein or Hamming distance
/// 
/// @param seq A string vector of sequences.
/// @param count An integer vector corresponding to sequence counts.
/// @param count_threshold A Integer Threshold to consider a sequence as true barcode.
/// @param dist_threshold A Integer Distance threshold for clustering.
/// @param depth_fold_threshold Fold change threshold for merging.
/// @param dist_method Distance method: 2 for Levenshtein, otherwise Hamming.
/// @param cost A vector of three integers for Levenshtein distance costs, the order is insert, delete, replace.
/// 
/// @return A list with two data frames:
/// - `seq_freq_tab`: Corrected sequence frequency table.
/// - `link_tab`: Clustering process record.
#[extendr]
pub fn seq_correct(
    seq: Vec<String>, 
    count: Vec<i32>, 
    count_threshold: i32, 
    dist_threshold: i32, 
    depth_fold_threshold: f64, 
    dist_method: i32, 
    cost: Vec<i32> // insert, delete, replace
) -> Robj {

    let mut lev = Levenshtein::new();
    lev.set_insert_cost(cost[0]);
    lev.set_delete_cost(cost[1]);
    lev.set_replace_cost(cost[2]);

    // Sort the sequences by count
    let mut old_seq: Vec<(String, i32)> = seq.into_iter().
        zip(count).collect();
    old_seq.sort_by(|a, b| b.1.cmp(&a.1));

    // Result Vector
    let mut res: Vec<(String, i32)> = Vec::new();

    // Clustering Tracking
    let mut remove_from: Vec<String> = Vec::new();
    let mut remove_from_size: Vec<i32> = Vec::new();
    let mut remove_by: Vec<String> = Vec::new();
    let mut remove_by_size: Vec<i32> = Vec::new();

    while !old_seq.is_empty() {
        if old_seq.len() == 1 {
            res.push(old_seq.pop().unwrap());
            break;
        }

        let tiptoe = old_seq.last().unwrap();

        if tiptoe.1 >= count_threshold {
            res.append(&mut old_seq);
            break;
        } else {
            let mut min_dist = i32::MAX;
            let mut min_it = 0;

            for (idx, branch) in old_seq.iter().enumerate().take(old_seq.len() - 1) {
                let h_dist = if dist_method == 2 {
                    crates_levenshtein_distance(
                        &branch.0,
                        &tiptoe.0,
                        &lev
                    )
                } else {
                    hamm_dist(&branch.0, &tiptoe.0)
                };

                if h_dist < min_dist {
                    min_dist = h_dist;
                    min_it = idx;
                }

                if h_dist == 1i32 {
                    break;
                }
            }

            let branch = &old_seq[min_it];

            if min_dist <= dist_threshold && (branch.1 as f64 / tiptoe.1 as f64) >= depth_fold_threshold {
                remove_from.push(tiptoe.0.clone());
                remove_from_size.push(tiptoe.1);
                remove_by.push(branch.0.clone());
                remove_by_size.push(branch.1);
                old_seq.pop();
            } else {
                res.push(old_seq.pop().unwrap());
            }
        }
    }

    // Sequence Frequency Table
    let res_seq: Vec<String> = res.iter().map(|x| x.0.clone()).collect();
    let res_count: Vec<i32> = res.iter().map(|x| x.1).collect();
    let seq_freq_tab = data_frame!(
        barcode_seq = res_seq,
        count = res_count
    );

    // Link Table
    let link_tab = data_frame!(
        seq_from = remove_from,
        seq_to = remove_by,
        from_size = remove_from_size,
        to_size = remove_by_size
    );

    list!(seq_freq = seq_freq_tab, link_table = link_tab).into()
}

extendr_module! {
    mod lib_clustering;
    fn seq_correct;
}
