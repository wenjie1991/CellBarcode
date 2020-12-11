#include <Rcpp.h>
using namespace Rcpp;

//' Hammer Distance
//' 
//' This function return hammer distance between two string.
//' If the two string do not have the same length, it will
//' return 999.
//'
//' @param s1, s2 two string
int hammer_dist(std::string s1, std::string s2) {
    if (s1.length() != s2.length()) {
        return 999;
    }
    int res = 0;
    for (int i=0; i < s1.length(); i++) {
        if (s1[i] != s2[i]) {
            res++;
        }
    }
    return res;
}

bool sortbycount(const std::pair<std::string, int> &a, const std::pair<std::string, int> &b) {
    return a.second > b.second;
}
    


//' Correct depth by UMI
//' 
//' This function will merge the UMIs by using the 
//' hammer distance. If two UMIs have hammer distance
//' no more than 1, only the UMI with more reads
//' will be kept.
//' 
//' This function will return the corrected UMI list.
//'
//' @param umi An String vector with UMI
//' @param count An integer vector with the same order and 
//' length of UMI
//' @export
// [[Rcpp::export]]
DataFrame seq_correct(std::vector<std::string> seq, IntegerVector count, int count_threshold, int hammer_dist_threshold) {
    
    // tree plot without fixed, all the nodes are candidates
    std::vector<std::pair<std::string, int>> cand;
    // a object to store the corrected sequence
    std::vector<std::pair<std::string, int>> res;
    
    for (auto i=0; i<seq.size(); i++) {
        cand.push_back(std::make_pair(seq[i], count[i]));
    }

    std::sort(cand.begin(), cand.end(), sortbycount);

    bool flag_all_meet_threshold = false;

    // while the flag does not true, which means remain some barcode
    // with count <  count_threshold
    while (!flag_all_meet_threshold & (!cand.empty())) {
        // tiptoe, the barcode with smallest count
        std::vector<std::pair<std::string, int>>::iterator tiptoe = cand.end() - 1;

        //std::cout<<tiptoe->first<<std::endl;

        // if the tiptoe pass the threshold stop the program 
        if (tiptoe->second >= count_threshold) {
            flag_all_meet_threshold = true;

            // store the candidates sequence
            res.insert(res.begin(), cand.begin(), cand.end());

            // stop cutting tiptoe 
            break;

        // else compare the tiptoe to the branch 
        } else {

            bool flag_is_connetcted = false;
        
            // if only one nodes left in the candidates list, stop
            if (cand.size() == 1) {
                res.insert(res.begin(), cand.begin(), cand.end());
                break;
            } 

            // find out if tiptoe connect to branch
            std::vector<std::pair<std::string, int>>::iterator it = cand.begin();
            while (it != tiptoe) {
                int h_dist = hammer_dist(it->first, tiptoe->first);
                // if the tiptoe connect to branch
                if (h_dist <= hammer_dist_threshold) {
                    // add the tiptoe to the branch node
                    it->second += tiptoe->second;
                    //std::cout<<it->second<<std::endl;
                    // remove the tiptoe
                    cand.pop_back();
                    // flag: if tiptoe is connect to branch
                    flag_is_connetcted = true;

                    // if only one nodes left in the candidates list, stop
                    if (!cand.empty()) {
                    //if (false) {

                        // update the order of the nodes by updated count
                        std::vector<std::pair<std::string, int>>::iterator it_in_rank = it;
                        while (it_in_rank != cand.begin()) {
                            it_in_rank--;
                            if (it_in_rank->second >= it->second) {
                                break;
                            }
                        }

                        if (it != it_in_rank) {
                            cand.insert(it_in_rank + 1, *it);
                            cand.erase(it);
                        }
                    } 

                    break;
                }
                it++;
            }

            // if tiptoe is not connect to any branch
            if (!flag_is_connetcted) {
                res.push_back(*tiptoe);
                cand.pop_back();
            }

        }
    }

    std::vector<std::string> res_seq;
    std::vector<int> res_count;

    for (auto i=0; i<res.size(); i++) {
        res_seq.push_back(res[i].first);
        res_count.push_back(res[i].second);
    }
    
    DataFrame df = DataFrame::create(Named("barcode_seq") = res_seq, Named("count") = res_count);
    return df;
}

// while (umi.size() != 0) {
//     int i = which_max(count);
//     std::string umi_current = umi[i];
//     res.push_back(umi_current);
//     
//     umi.erase(umi.begin() + i);
//     count.erase(count.begin() + i);
//     
//     int size = umi.size();
//     for (int i=0; i<size; i++) {
//         int h_distance = hammer_dist(umi_current, umi[i]);
//         if (h_distance == 1) {
//             umi.erase(umi.begin() + i);
//             count.erase(count.begin() + i);
//         }
//     } 
// }


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
