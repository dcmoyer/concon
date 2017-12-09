
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

