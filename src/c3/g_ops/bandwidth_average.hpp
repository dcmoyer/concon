
#ifndef C3_BANDWIDTH_AVE_HPP
#define C3_BANDWIDTH_AVE_HPP

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include "c3/csv.h"
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

