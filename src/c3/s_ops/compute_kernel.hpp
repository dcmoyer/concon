//Copyright 2017 Daniel Moyer
//
//Permission is hereby granted, free of charge, to any person obtaining a
//copy of this software and associated documentation files (the "Software"),
//to deal in the Software without restriction, including without limitation
//the rights to use, copy, modify, merge, publish, distribute, sublicense,
//and/or sell copies of the Software, and to permit persons to whom the
//Software is furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in
//all copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
//OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
//ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
//OTHER DEALINGS IN THE SOFTWARE.
//


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
