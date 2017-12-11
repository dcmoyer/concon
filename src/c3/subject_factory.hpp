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


#ifndef C3_SUBJ_FACTORY_HPP
#define C3_SUBJ_FACTORY_HPP


#include "c3/subject.hpp"
#include "c3/s_ops/integrate_label.hpp"
#include "c3/s_ops/bandwidth_estimation.hpp"
#include "c3/s_ops/compute_kernel.hpp"
#include "c3/s_ops/effective_connectivity.hpp"
#include "c3/s_ops/surface_intersect.hpp"

#include <memory>
#include <string>
#include <iostream>


namespace c3 {

class Subject_Factory {
public:
  std::shared_ptr<c3::Subject> make_subject(std::string type);
};



}//end of namespace c3

#endif

