
#ifndef C3_CREATE_GROUP_MASK_HPP
#define C3_CREATE_GROUP_MASK_HPP

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <omp.h>

namespace c3 {

  int create_group_mask( std::vector<std::string> filenames,
    std::string output, bool verbose=false);

}//end of namespace c3

#endif

