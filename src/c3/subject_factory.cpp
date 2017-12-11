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


#include "c3/subject_factory.hpp"

namespace c3{

std::shared_ptr<c3::Subject> Subject_Factory::make_subject(std::string type){

  if(type == "Bandwidth_Estimation"){
    return std::shared_ptr<c3::Subject>(new c3::Bandwidth_Estimation());
  }

  if(type == "Compute_Kernel"){
    return std::shared_ptr<c3::Subject>(new c3::Compute_Kernel());
  }

  if(type == "Effective_Connectivity"){
    return std::shared_ptr<c3::Subject>(new c3::Effective_Connectivity());
  }

  if(type == "Integrate_Label"){
    return std::shared_ptr<c3::Subject>(new c3::Integrate_Label());
  }

  if(type == "Surface_Intersect"){
    return std::shared_ptr<c3::Subject>(new c3::Surface_Intersect());
  }


  std::cout << "ERROR: That type has not yet been implemented. Returns NULL"
    << std::endl;
  return NULL;

}

}//end of namespace c3




