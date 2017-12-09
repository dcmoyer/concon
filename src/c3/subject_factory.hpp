
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

