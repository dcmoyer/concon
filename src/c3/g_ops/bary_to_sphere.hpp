
#ifndef C3_BARY_TO_SPHERE_HPP
#define C3_BARY_TO_SPHERE_HPP

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <utility>
#include <stdio.h>
#include <unordered_map>

#include "MeshLib/Solid.h"
#include "MeshLib/iterators.h"
#include "c3/subject.hpp"

namespace c3 {

  int bary_to_sphere(
    std::unordered_map<std::string,std::string> params,
    std::string subj=""
  );

}//end of namespace c3

#endif


