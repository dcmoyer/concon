
#ifndef C3_MERGE_MARGINALS_HPP
#define C3_MERGE_MARGINALS_HPP

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>

namespace c3 {

  int merge_marginals( std::vector<std::string> subj_list,
    std::unordered_map<std::string,std::string> params);

}//end of namespace c3

#endif


