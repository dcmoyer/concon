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


#ifndef C3_BANDWIDTH_AVE_HPP
#define C3_BANDWIDTH_AVE_HPP

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
//#include "c3/csv.h"
#include "c3/subject.hpp"

namespace c3 {

int Bandwidth_Average( std::vector<std::string> subj_list,
  std::unordered_map<std::string,std::string> params);

//
//  Bandwidth Average
//

class Bandwidth_Average_Obj {
public:
  Bandwidth_Average_Obj(std::vector<std::string> files, int rs,
    double r_min, double r_max, std::string out);
  ~Bandwidth_Average_Obj();

  int load_files();
  int save_to_csv();

  //uses mean
  int compute_optimum();

private:
  std::vector<std::string> filenames;
  std::string output_name;
  int range_size;
  double range_min;
  double range_max;
  double min_obj;
  int min_obj_idx;
  double** sigma_objective;

};

}//end of namespace c3

#endif

