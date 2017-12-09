
#ifndef C3_CREATE_GROUP_LABEL_MASK_HPP
#define C3_CREATE_GROUP_LABEL_MASK_HPP

#include <vector>
#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <sstream>
#include "MeshLib/Solid.h"
#include "MeshLib/iterators.h"

namespace c3 {

  int create_group_label_mask( std::vector<std::string> subj_list,
    std::unordered_map<std::string,std::string> params);

}//end of namespace c3

#endif

