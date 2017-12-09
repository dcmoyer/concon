
///
/// :::Integrate Kernels:::
///
/// This class integrates kernels for a specific label or all (-1)
///

#include "c3/subject.hpp"

#ifndef C3_INTEGRATE_LABEL_HPP
#define C3_INTEGRATE_LABEL_HPP

namespace c3 { 

class Integrate_Label : public c3::Subject {
public:
  Integrate_Label();
  ~Integrate_Label();

  int set(std::unordered_map<std::string, std::string> params);

  int subj_specific_load(std::string subj);

  bool check_req();  

  int set_display_mesh(std::string lh_file, std::string rh_file);

  int save_file(std::string filename);

  int load_file(std::string filename);

  int run();

protected:
  //add new data here

  std::vector<int> target_regions;
  std::vector<int> lh_labels, rh_labels;

  std::string label_prefix;
  std::string label_lh_postfix;
  std::string label_rh_postfix;

  std::string subj_id;
  std::string out_prefix;
  std::string out_postfix;

  bool mean_flag;

};

}//end of namespace c3

#endif

