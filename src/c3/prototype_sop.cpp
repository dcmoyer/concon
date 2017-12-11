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


#include "c3/prototype_sop.hpp"

namespace c3{

int OP_NAME::set(std::unordered_map<std::string, std::string> params){

  //set vars here
  //
  for( std::pair<std::string,std::string> p : params ){

  }

}

int OP_NAME::subj_specific_load(std::string subj){
  return 0;
}

bool OP_NAME::check_req();
  
int OP_NAME::run(){

  if(!check_req()){
    //condition?
    return 1;
  }

  //do operation here.

}

int OP_NAME::save_file(std::string filename){

}

int OP_NAME::load_file(std::string filename){

}

}//end of namespace c3

