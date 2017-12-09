
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




