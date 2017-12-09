
#include "c3/g_ops/bandwidth_average.hpp"

namespace c3 {

int Bandwidth_Average(
  std::vector<std::string> subj_list,
  std::unordered_map<std::string,std::string> params
){

  std::vector<std::string> filenames;
  for(std::string & s : subj_list)
    filenames.push_back(params["SAVE_Bandwidth_Estimation_prefix"] +
      s +
      params["SAVE_Bandwidth_Estimation_postfix"]);

  int rs = -1;
  double r_min = -1,r_max = -1;
  for(std::pair<std::string,std::string> p : params){

    if(p.first == "range_size"){
      rs = c3::string_to_int(p.second);
    }
    if(p.first == "range_min"){
      r_min = c3::string_to_double(p.second);
    }
    if(p.first == "range_max"){
      r_max = c3::string_to_double(p.second);
    }
 
  }

  if( rs == -1 ){
    std::cerr << "Error: Range Size Not Set, Aborting." << std::endl;
    return 1;
  }

  std::string outfile = params["LOAD_sigma_file"];

  c3::Bandwidth_Average_Obj BA(filenames,rs,r_min,r_max,outfile);

  BA.load_files();
  BA.compute_optimum();
  BA.save_to_csv();

  return 0;
}

Bandwidth_Average_Obj::Bandwidth_Average_Obj(std::vector<std::string> files,
  int rs, double r_min, double r_max, std::string out){
  filenames = files;
  range_size = rs;
  range_min = r_min;
  range_max = r_max;
  output_name = out;

  int s = filenames.size();
  sigma_objective = new double*[s];
  for(int i = 0; i < s; ++i){
    sigma_objective[i] = new double[rs];
  }

  min_obj_idx = -1;
  min_obj = -1;
}

Bandwidth_Average_Obj::~Bandwidth_Average_Obj(){

  if(sigma_objective != NULL){
    int s = filenames.size();
    for(int i = 0; i < s; ++i)
      delete[] sigma_objective[i];
    delete[] sigma_objective;
  }

}

int Bandwidth_Average_Obj::load_files(){

  int s = filenames.size();
  for(int i = 0; i < s; ++i){
    std::ifstream input(filenames[i].c_str());
    for(int j = 0; j < range_size; ++j){
      input >> sigma_objective[i][j];
    }
    input.close();
  }

  return 0;
}

int Bandwidth_Average_Obj::compute_optimum(){
  double * sum_obj = new double[range_size];

  for(int i = 0; i < range_size; ++i){
    sum_obj[i] = 0;
  }

  int s = filenames.size();
  for(int i = 0; i < s; ++i){
    for(int j = 0; j < range_size; ++j){
      sum_obj[j] += sigma_objective[i][j];
    }
  }

  double running_min = sum_obj[0];
  int running_min_idx = -1;
  for(int i = 0; i < range_size; ++i){
    if(sum_obj[i] < running_min){
      running_min = sum_obj[i];
      running_min_idx = i;
    }
  }

  if(running_min_idx == -1){
    std::cerr << "ERROR: Optimization error, Aborting" << std::endl;
    return -1;
  }

  min_obj = running_min;
  min_obj_idx = running_min_idx;
  

  return running_min_idx;
}

int Bandwidth_Average_Obj::save_to_csv(){

  std::ofstream out(output_name.c_str());

  out << min_obj_idx * (range_max - range_min)/range_size + range_min;

  out.close();

  return 0;
}

}//end of namespace c3


