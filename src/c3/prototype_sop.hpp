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

