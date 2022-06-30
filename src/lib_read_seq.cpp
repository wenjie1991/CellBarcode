#include <Rcpp.h>
#include <algorithm>
#include <unordered_map>
#include <string>
#include <zlib.h>
using namespace Rcpp;

//' Read the Fastq file and output frequency table
//'
//' This function read in fastq.gz file and output the sequences frequency
//' data.frame with two columns, freq and seq.
// read the fastq.gz
// write the fastq.gz

// [[Rcpp::export]]
DataFrame read_fastq_gz(std::string in_file_path) {

    // input
    gzFile infile = gzopen(in_file_path.c_str(), "r");

    // output

    // buffer
    char line_title[0x10000];
    char line_seq[0x10000];
    char line_sep[0x10000];
    char line_quality[0x10000];
    std::unordered_map<std::string, long long> seq_map;

    while (!gzeof(infile)) {

        // n++;
        // if (n > 5) {
        //     break;
        // }

        gzgets(infile, line_title, 1000);
        gzgets(infile, line_seq, 1000);
        gzgets(infile, line_sep, 1000);
        gzgets(infile, line_quality, 1000);
        line_seq[strcspn(line_seq, "\n")] = 0;

        if (seq_map.find(line_seq) == seq_map.cend()) {
            seq_map.insert(std::pair<std::string, long long>(line_seq, 1));
        } else {
            seq_map[line_seq]++;
        }
    }

    gzclose(infile);

    IntegerVector freq(seq_map.size());
    CharacterVector seq(seq_map.size());

    auto i = 0;
    for (auto it = seq_map.begin(); it != seq_map.end(); it++) {
        seq[i] = it->first;
        freq[i] = it->second;
        i++;
    }

    return DataFrame::create( _["freq"] = freq , _["seq"] = seq );
}

// [[Rcpp::export]]
DataFrame read_fastq(std::string in_file_path) {

    // input
    FILE * infile = fopen(in_file_path.c_str(), "r");

    // output

    // buffer
    char line_title[0x10000];
    char line_seq[0x10000];
    char line_sep[0x10000];
    char line_quality[0x10000];
    std::unordered_map<std::string, long long> seq_map;

    while (!feof(infile)) {

        // n++;
        // if (n > 5) {
        //     break;
        // }

        fgets(line_title, 1000, infile);
        fgets(line_seq, 1000, infile);
        fgets(line_sep, 1000, infile);
        fgets(line_quality, 1000, infile);
        line_seq[strcspn(line_seq, "\n")] = 0;

        if (seq_map.find(line_seq) == seq_map.cend()) {
            seq_map.insert(std::pair<std::string, long long>(line_seq, 1));
        } else {
            seq_map[line_seq]++;
        }
    }

    fclose(infile);

    IntegerVector freq(seq_map.size());
    CharacterVector seq(seq_map.size());

    auto i = 0;
    for (auto it = seq_map.begin(); it != seq_map.end(); it++) {
        seq[i] = it->first;
        freq[i] = it->second;
        i++;
    }

    return DataFrame::create( _["freq"] = freq , _["seq"] = seq );
}
