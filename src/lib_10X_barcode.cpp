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

//' @param sam A string, define the un-mapped sequences 
//' @param pattern A string, define the regular expression to match the barcode
//' sequence. The barcode sequence should be in the first catch. Please see the
//' \code{\link[CellBarcode]{bc_extract}} for detail.
//' @param cell_barcode_tag A string, define the tag of 10X cell barcode field in sam
//' file. The default is "CR".
//' @param umi_tag A string, define the tag of UMI field in the sam file.
//' @return 
//' A data.frame with 4 columns:
//' \enumerate{
//'   \item \code{cell_barcode}: 10X cellular barcode.
//'   \item \code{umi}: UMI sequence.
//'   \item \code{barcode_seq}: lineage barcode.
//'   \item \code{count}: reads count.
// [[Rcpp::export]]
DataFrame parse_10x_scSeq(
    std::string in_file_path, 
    std::string regex_str, 
    std::string cell_barcode_tag = "CR",
    std::string umi_tag = "UR"
) {

    // INPUT
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
            if (parts[j].substr(0, 2) == cell_barcode_tag) {
                cell_barcode = parts[j].substr(5);
            } else if (parts[j].substr(0, 2) == umi_tag) {
                umi_seq = parts[j].substr(5);
            }
        }

        if (boost::regex_search(seq, sm, str_expr)) {

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
            _["barcode_seq"] = barcode_v,
            _["count"] = count_v
        );
}
