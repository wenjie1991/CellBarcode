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


//' Sequence clustering
//' 
//' This function will merge the UMIs by using the 
//' hammer distance. If two UMIs have hammer distance
//' no more than 1, only the UMI with more reads
//' will be kept.
//' 
//' This function will return the corrected UMI list.
//'
//' @param seq A string vector
//' @param count An integer vector with the same order and length of UMI
//' @param count_threshold An integer the maximum barcodes number
//' @param hammer_dist_threshold An integer the hammer distance threshold
// [[Rcpp::export]]
List seq_correct(std::vector<std::string> seq, IntegerVector count, int count_threshold, int hammer_dist_threshold) {

  // candidates: all the nodes are candidates
  std::vector<std::pair<std::string, int>> cand;

  // results: a object to store the corrected sequence
  std::vector<std::pair<std::string, int>> res;


  // the small node be merged
  std::vector<std::string> merge_from;
  std::vector<int> merge_from_size;
  // the big node be merged
  std::vector<std::string> merge_to;
  std::vector<int> merge_to_size;

  // Sort by the frequency of seq
  for (auto i=0; i<seq.size(); i++) {
    cand.push_back(std::make_pair(seq[i], count[i]));
  }
  std::sort(cand.begin(), cand.end(), sortbycount);

  //bool flag_all_meet_threshold = false;

  // while the flag does not true, which means remain some barcode
  // with count <  count_threshold
  while (!cand.empty()) {

    // if only one nodes left in the candidates list, stop
    if (cand.size() == 1) {
      res.insert(res.begin(), cand.begin(), cand.end());
      break;
    } 

    // tiptoe: the barcode with smallest count
    std::vector<std::pair<std::string, int>>::iterator tiptoe = cand.end() - 1;

    //std::cout<<tiptoe->first<<std::endl;

    // if the tiptoe pass the threshold stop the program 
    if (tiptoe->second >= count_threshold) {
      // flag_all_meet_threshold = true;

      // store the candidates sequence
      res.insert(res.begin(), cand.begin(), cand.end());

      // stop cutting tiptoe 
      break;

      // else compare the tiptoe to the branch 
    } else {

      // is tiptoe is connect to big nodes
      bool flag_is_connetcted = false;
      int min_dist = 2147483646;
      int h_dist;
      std::vector<std::pair<std::string, int>>::iterator min_it;


      // find out if tiptoe connect to branch
      std::vector<std::pair<std::string, int>>::iterator it = cand.begin();
      while (it != tiptoe) {
        h_dist = hammer_dist(it->first, tiptoe->first);
        // record the min_dist
        if (h_dist < min_dist) {
          min_dist = h_dist;
          min_it = it;
        }
        // if the branch is the nearest one to tiptoe 
        if (h_dist == 1) {
          break;
        } 
        it++;
      }

      if (min_dist <= hammer_dist_threshold) {
        // add the tiptoe to the branch node
        min_it->second += tiptoe->second;
        //std::cout<<min_it->second<<std::endl;
        // record the nodes
        merge_from.push_back(tiptoe->first);
        merge_from_size.push_back(tiptoe->second);
        merge_to.push_back(min_it->first);
        merge_to_size.push_back(min_it->second);

        // remove the tiptoe
        cand.pop_back();
        // flag: if tiptoe is connect to branch
        flag_is_connetcted = true;

        // update the order of the nodes using updated count
        std::vector<std::pair<std::string, int>>::iterator it_in_rank = min_it;
        while (it_in_rank != cand.begin()) {
          it_in_rank--;
          if (it_in_rank->second >= min_it->second) {
            break;
          }
        }

        if (min_it != it_in_rank) {
          cand.insert(it_in_rank + 1, *it);
          cand.erase(min_it);
        }
      } else {
        // if tiptoe is not connect to any branch
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

  // seqeunce frequency table
  DataFrame seq_freq_tab = DataFrame::create(Named("barcode_seq") = res_seq, Named("count") = res_count);

  // link table
  DataFrame link_tab = DataFrame::create(Named("seq_from") = merge_from, Named("seq_to") = merge_to, _["from_size"] = merge_from_size, _["to_size"] = merge_to_size);

  List L = List::create(Named("seq_freq") = seq_freq_tab, _["link_table"] = link_tab);

  return L;
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
