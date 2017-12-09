
#include "c3/s_ops/integrate_label.hpp"

namespace c3{

Integrate_Label::Integrate_Label() : c3::Subject() {
  r_kernel = true;
  mean_flag = false;
  return;
}

Integrate_Label::~Integrate_Label(){
  return;
}

int Integrate_Label::set(
  std::unordered_map<std::string, std::string> params){

  //set vars here
  //
  for( std::pair<std::string,std::string> p : params ){
    if(p.first == "verbose"){
      verbose = true;
    }

    if(p.first == "regions"){
      std::string temp = p.second;
      std::stringstream ss;
      ss << temp;
      while(!ss.eof()){
        int temp_target;
        ss >> temp_target;
        target_regions.push_back(temp_target);
        ss.get();

        if(verbose)
          std::cout << temp_target << " ";
      }

      if(verbose)
        std::cout << std::endl;
    }

    if(p.first == "label_prefix"){
      label_prefix = p.second;
    }

    if(p.first == "label_lh_postfix"){
      label_lh_postfix = p.second;
    }

    if(p.first == "label_lh_postfix"){
      label_rh_postfix = p.second;
    }

    if(p.first == "out_prefix"){
      out_prefix = p.second;
    }

    if(p.first == "out_postfix"){
      out_postfix = p.second;
    }

    if(p.first == "mean_flag"){
      mean_flag = true;
    }
    //end of param reading
  }
  return 0;
}

int Integrate_Label::subj_specific_load(std::string subj){

  subj_id = subj;

  if(verbose)
    std::cout << "Loading labels" << std::endl;

  std::string lh_file = label_prefix + subj + label_lh_postfix;
  std::ifstream input_lh_file(lh_file.c_str());
  int temp;
  int max = -1;
  std::string rh_file = label_prefix + subj + label_rh_postfix;
  std::ifstream input_rh_file(rh_file.c_str());

  while(input_rh_file >> temp){
    if(temp > max)
      max = temp;
    rh_labels.push_back(temp);
  }
  input_rh_file.close();

  while(input_lh_file >> temp){
    if(temp >= 0){
      lh_labels.push_back(temp + max);
    } else { 
      lh_labels.push_back(temp);
    }
  }
  input_lh_file.close();


  if(verbose){
    std::cout << "rh label size " << rh_labels.size() << std::endl;
    std::cout << "lh label size " << lh_labels.size() << std::endl;
  }

  if((int) (rh_labels.size() + lh_labels.size()) != kernel->size()){
    std::cout << "Error: kernel does not match label sizes, aborting."
      << std::endl;
    exit(1);
  }

  return 0;
}

bool Integrate_Label::check_req(){
  return true;
}
  
int Integrate_Label::run(){

  if(!check_req()){
    //condition?
    return 1;
  }

  if(verbose)
    std::cout << "Kernel Size " << kernel->size() << std::endl;

  std::vector<int> target_idxs;
  if(target_regions[0] == -1919){
    //if signal given to use all labels
    target_regions.clear();

    //start with rh
    target_regions = rh_labels;
    std::unique(target_regions.begin(),target_regions.end());

    //copy over lh
    std::vector<int> second_unique = lh_labels;
      std::unique(second_unique.begin(),second_unique.end());
    target_regions.reserve(target_regions.size() + second_unique.size());
    target_regions.insert(target_regions.end(),
      second_unique.begin(), second_unique.end());
  }

  //marginalize like we're Western Society
  int temp;
  for(int& region : target_regions){
    if(region == -2){
      //full_marginal_conn

      std::ofstream out(out_prefix + subj_id + "_full"  + out_postfix);
      int s_rh = rh_labels.size();
      int s_lh = lh_labels.size();

      if(!mean_flag){
        for(int i = 0; i < s_rh + s_lh; ++i){
          out << kernel->row_sum(i) << std::endl;
        }
      } else {
        for(int i = 0; i < s_rh + s_lh; ++i){
          temp = kernel->row_num_nonzeros(i);
          if(temp > 0){
            out << kernel->row_sum(i)/((double) temp) << std::endl;
          } else {
            out << 0.0 << std::endl;
          }
        }
      }

      out.close();
    } else {

      // target one region only

      std::ofstream out(out_prefix + subj_id + "_" +
        std::to_string(region) + out_postfix);

      target_idxs.clear();
      int s_rh = rh_labels.size();
      int s_lh = lh_labels.size();
      for( int i = 0; i < s_rh; ++i ){
        if(rh_labels[i] == region)
          target_idxs.push_back(i);
      }
      if(target_idxs.size() == 0){
        for( int i = 0; i < s_lh; ++i ){
          if(lh_labels[i] == region)
            target_idxs.push_back(i + s_rh);
        }
      }
      for(int i = 0; i < s_rh + s_lh; ++i){
        out << kernel->ordered_array_sum(i, target_idxs) << std::endl;
      }
      out.close();

    }
  }

  return 0;
}

int Integrate_Label::save_file(std::string filename){
  return -1;
}

int Integrate_Label::load_file(std::string filename){
  return -1;
}

}//end of namespace c3

