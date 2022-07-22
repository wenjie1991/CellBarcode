#include <cstring>
#include <string>
#include <unordered_map>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

struct Barcode {
    std::string cell_barcode;
    std::string UMI;
    std::string barcode;
    long count;
};

// [[Rcpp::export]]
DataFrame parse_10x_bam(std::string in_file_path, std::string regex_str) {

    // INPUT
    // std::string in_file_path = "./data/test_sub.sam";
    // std::string in_file_path = "./data/test.sam";
    // std::regex str_expr ("CGTCAACTAGAACACTCGAGATCAG(.*)TGTGGTATGATGTATCATCTGAGTA");
    // std::regex str_expr (".*AGATCAG(.*)TGTGGTAT.*");
    // boost::regex str_expr (".*AGATCAG(.*)TGTGGTAT.*");
    boost::regex str_expr (regex_str);
    
    FILE * infile = fopen(in_file_path.c_str(), "r");

    char line[0x10000];
    char seq[0x10000];
    std::string cell_barcode;
    std::string umi_seq;
    std::string barcode;

    // std::smatch sm;
    boost::cmatch sm;

    std::unordered_map<std::string, Barcode> seq_map;

    int i = 0;
    std::string key;

    while (!feof(infile)) {
        fgets(line, 1000, infile);
        line[strcspn(line, "\n")] = 0;
        i++;
        // if (i == 1000000) {
            // break;
        // }
        std::vector<std::string> parts;
        split(parts, line, boost::is_any_of("\t"));
        // seq = parts[9];
        strcpy(seq, parts[9].c_str());
        for (int j=11; j<parts.size(); j++) {
            if (parts[j].substr(0, 2) == "CR") {
                cell_barcode = parts[j].substr(5);
            } else if (parts[j].substr(0, 2) == "UR") {
                umi_seq = parts[j].substr(5);
            }
        }

        if (boost::regex_match(seq, sm, str_expr)) {

            // std::cout 
            //     << "seq:" << seq << std::endl
            //     << "cell_barcode:" << cell_barcode << std::endl
            //     << "umi:" << umi_seq << std::endl
            //     << "match:" << sm[1] << std::endl
            //     << std::endl;

            barcode = sm[1];

            key = cell_barcode + umi_seq + barcode;

            auto iter = seq_map.find(key);

            if (iter != seq_map.end()) {
                iter->second.count++;
            } else {
                Barcode barcode_st;
                barcode_st.cell_barcode = cell_barcode;
                barcode_st.UMI = umi_seq;
                barcode_st.barcode = barcode;
                barcode_st.count = 1;

                seq_map.insert(std::make_pair(key, barcode_st));
            }
        }
        // std::cout<<line;
    }

    fclose(infile);

    CharacterVector barcode_v (seq_map.size());
    CharacterVector cell_barcode_v (seq_map.size());
    CharacterVector umi_v (seq_map.size());
    IntegerVector count_v (seq_map.size());

    int j = 0;
    for (auto iter = seq_map.begin(); iter != seq_map.end(); iter++) {
        barcode_v[j] = iter->second.barcode;
        cell_barcode_v[j] = iter->second.cell_barcode;
        umi_v[j] = iter->second.UMI;
        count_v[j] = iter->second.count;

        j++;
    }
    
    return DataFrame::create(
            _["cell_barcode"] = cell_barcode_v,
            _["umi"] = umi_v,
            _["count"] = count_v
        );
}
