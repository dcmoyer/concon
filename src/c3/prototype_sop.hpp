
///
/// :::PROTOTYPE! DO NOT ACTUALLY USE:::
///
/// This is just an example.
///

#include <string>
#include <iostream>
#include <vector>
#include <utility>
#include <unordered_map>

#include "c3/subject.hpp"

#ifndef C3_OP_NAME_HPP
#define C3_OP_NAME_HPP

namespace c3 { 

class OP_NAME : public c3::Subject {
public:
  OP_NAME() : c3::Subject(){;} //usually blank?
  ~OP_NAME(){;} //usually blank as well? Unless you defined vars.

  int set(std::unordered_map<std::string, std::string> params);

  int subj_specific_load(std::string subj);

  bool check_req();  

  int save_file(std::string filename);

  int load_file(std::string filename);

  int run();

protected:
  //add new data here

};

}//end of namespace c3

#endif

