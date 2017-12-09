
//
//  Compute_Kernel
//
//  This class encapsulates the actual computation of the kernel.
//
//Notes about the output format
//  -RIGHT comes first. Because Boris said so like 20 years ago.

#ifndef C3_COMPUTE_KERNEL_HPP
#define C3_COMPUTE_KERNEL_HPP

#include <iomanip>
#include "c3/subject.hpp"
#include "c3/sigma_opt.hpp"
#include "utils/quad_tree.hpp"
#include "MeshLib/Solid.h"

namespace c3 { 

class Compute_Kernel : public c3::Subject {
public:
  Compute_Kernel();
  ~Compute_Kernel();

  int set(std::unordered_map<std::string, std::string> params);

  int subj_specific_load(std::string subj);

  bool check_req();

  int save_file(std::string filename);

  int load_file(std::string filename);
 
  int run();

protected:
  int N_output;
  double ** output;
  double epsilon;
  double final_thold;
  double *temp_harmonics_0;
  double *temp_harmonics_1;

  std::string trk_weight_prefix;
  std::string trk_weight_postfix;
  bool weights_set;
  double *trk_weights;

};

}//end of namespace c3




#endif
