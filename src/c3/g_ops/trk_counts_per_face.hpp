
#ifndef C3_G_OPS_TRK_ENDPOINTS_PER_FACE
#define C3_G_OPS_TRK_ENDPOINTS_PER_FACE

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include <unordered_map>
#include <stdio.h>
#include "MeshLib/Solid.h"
#include "MeshLib/iterators.h"
#include "c3/subject.hpp"

namespace c3 {

int trk_endpoints_per_face(
  std::unordered_map< std::string, std::string > params
);

}//end of namespace c3



#endif



