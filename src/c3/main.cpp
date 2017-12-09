

#include <iostream>
#include <string>
#include <memory>
#include <unordered_map>
#include "utils/interface.hpp"
#include "c3/subject.hpp"
#include "c3/subject_factory.hpp"
#include "c3/s_ops/bandwidth_estimation.hpp"
#include "c3/g_ops/bandwidth_average.hpp"
#include "c3/g_ops/create_group_mask.hpp"
#include "c3/g_ops/merge_marginals.hpp"
#include "c3/g_ops/create_group_label_mask.hpp"
#include "c3/g_ops/trk_counts_per_face.hpp"
#include "c3/g_ops/bary_to_sphere.hpp"

void print_help(){
  std::cerr << "c3 main function: " << std::endl
    << "\targ1:\t[command]" << std::endl
    << "\targ2:\t--env-file" << std::endl
    << "\targ3:\t[env file.tsv]" << std::endl
    << std::endl

    << "Overrides:" << std::endl
    << "\targN:\t--VAR-TO-OVERRIDE" << std::endl
    << "\targ(N+1):\tVALUE" << std::endl
    << std::endl

    << "Notes:" << std::endl
    << "\tVAR chars '-' and '_' are interchangable besides the initial '--'."
    << std::endl
    << "\tCurrently each arg needs a value (flag style support to come)."
    << std::endl << std::endl
    << "\tSubjects can either be given using '--subj' or in a list with "
    << "'--subj-list'." << std::endl
    << "\tGrid handling can be passed using '--job-idx' or '--partition-num'."
    << std::endl
    << "\tZero Indexing is default, One Indexing can be used with "
    << "'--one-indexing'." << std::endl
    << "\tSetting '--max-partitions' will allow for different numbers of "
    << "jobs and subjects." << std::endl
    << "\tCould get weird tho." << std::endl
    << std::endl

    << "No args brings up a fun few paragraphs." << std::endl
    << std::endl;
}

void verbose_help(){

  std::cout << std::endl
    << "This is the main c3 binary." << std::endl << std::endl

    << "The '-h' flag brings up the actual help." << std::endl << std::endl

    << "Given a set of tracts (.trk) and a surface (.m) in the same space(!)"
    << std::endl
    <<  "it produces the tract intersections with that surface "
    "(as bary-centric coordinates)."
    << std::endl << std::endl

    << "Given a registered surface to a sphere, and tract intersection"
    "coordinates" << std::endl
    << "it produces the 'spherical gaussian' KDE estimate of the"
    << "connectivity function."
    << std::endl
    << "See Moyer et al. 2017 Medical Image Analysis for detals."
    << std::endl << std::endl

    << "It also does other stuff that I'm too lazy to document right now."
    << std::endl << std::endl

    << "It comes with ABSOLUTELY NO WARRENTY, and the author is not"
    << " responsible for your copy"
    << std::endl
    << "Or how you choose to use or modify it. Goodluck."
    << std::endl << std::endl
    << "D Moyer 171128"
    << std::endl << std::endl;

  return;
}

void print_version(){

  return;
}


//
//  Main function
//

int main(int argc, char *argv[]){

  if(argc < 2){
    verbose_help();
    return 1;
  }

  int verbose = 0;

  std::string command = std::string(argv[1]);

  std::unordered_map<std::string,std::string> params;
  std::string env_file;

  //
  // Load Parameters
  //

  //TODO: get params/settings from commands themselves.
  if(argc < 3){
    interface::load_params(argc, argv, params, env_file, print_help,
      print_version, 1);
  } else {
    interface::load_params(argc, argv, params, env_file, print_help,
      print_version, 2);
  }

  interface::load_env(params,env_file);

  if(params.find("verbose") != params.end()){
    verbose = 1;
  }

  //
  //  Subject Handling
  //

  //TODO: LOAD SUBJ FILE

  std::vector<std::string> subjects;
  if(params.find("subj") != params.end()) {
    //if subj given

    if(verbose > 0)
      std::cout << "[c3/main] " <<
        "Using subj " << params["subj"] << std::endl;

    subjects.push_back(params["subj"]);
  } else if (params.find("subj_list") != params.end()) {
    // if subj_file given 

    if(verbose > 0)
      std::cout << "[c3/main] " << 
        "Loading subj list from " << params["subj_file"] << std::endl;

    std::ifstream subj_input(params["subj_file"].c_str());
    std::string temp;
    while( subj_input >> temp ){
      subjects.push_back(temp);
    }
  } else {

    //otherwise if no subj nor subj_list
    //if not silent 
    if(params.find("silent") == params.end()){
      std::cerr << "[c3/main] " <<
        "WARNING: Neither subj nor subj_list set. " <<
        "Using empty string as subj. " << std::endl;
      if(verbose <= 0){
        std::cerr << 
          "(use --silent to silence this alert)." << std::endl;
      }
    }

    subjects.push_back("");
  }

  int num_subjects = subjects.size();

  //
  //  subj-wise grid compute handlers
  //    Note that other functions are still able to use grid job-arrays,
  //    they'll simply have to code it in themselves. The grid equivalent of
  //    SGE_TASK_ID needs to be passed via the --job_idx or --partition_num
  //    flags (they are equivalent).

  int partition_num;
  int max_partitions = num_subjects;

  if (params.find("job_idx") != params.end()) {
    partition_num = c3::string_to_int(params["job_idx"]);
  } else if (params.find("partition_num") != params.end()) {
    partition_num = c3::string_to_int(params["partition_num"]);
  } else {

    //just runs all subjs in subj_list, which could be just one, or even
    //the empty string subj

    partition_num = 0;
    max_partitions = 1;
  }

  if ((params.find("partition_num") != params.end() ||
    params.find("job_idx") != params.end()) &&
    params.find("max_partitions") != params.end()) {
    max_partitions = c3::string_to_int(params["max_partitions"]);
  }

  int one_indexing = 0;
  if(params.find("one_indexing") != params.end()){
    one_indexing = 1;
  }

  int part_size = 0;
  if(num_subjects % max_partitions == 0) {
    part_size = num_subjects / max_partitions;
  } else {
    part_size = (num_subjects / max_partitions) + 1;
  }

  std::vector<std::string>::const_iterator start_itr = subjects.begin()
    + (part_size * (partition_num-one_indexing));
  std::vector<std::string>::const_iterator end_itr = subjects.begin()
    + std::min(part_size * (partition_num - one_indexing + 1), num_subjects);

  subjects = std::vector<std::string>(start_itr, end_itr);

  //
  //  ops
  //

  //TODO:Deal with gopts
  if(command == "Bandwidth_Average"){

    if(verbose > 0)
      std::cout << "[c3/main]: " <<
        "Bandwidth Average" << std::endl;

    return c3::Bandwidth_Average( subjects, params );

  } else if (command == "create_group_mask"){
 
    std::vector<std::string> full_filenames;
    int s = subjects.size();
    for(int i = 0; i < s; ++i)
      full_filenames.push_back(params["SAVE_Compute_Kernel_prefix"]
        + subjects[i] + params["SAVE_Compute_Kernel_postfix"]);

    if(c3::create_group_mask( full_filenames,
      params["LOAD_mask_file"], true) != 0){
      std::cout << "Error in mask call, aborting" << std::endl;
      return 1;
    }

    return 0;
  } else if (command == "merge_marginals"){

    if(verbose > 0)
      std::cout << "Merging Marginals" << std::endl;

    c3::merge_marginals(subjects,params);

    return 0;
  } else if (command == "create_group_label_mask"){

    if(verbose > 0)
      std::cout << "Create Group Label Mask" << std::endl;

    c3::create_group_label_mask(subjects,params);

    return 0;
  } else if (command == "trk_counts_per_face") {

    if(verbose > 0)
      std::cout << "Trk Counts Per Face" << std::endl;

    c3::trk_endpoints_per_face( params );

    return 0;
  } else if (command == "bary_to_sphere") {

    if(verbose > 0)
      std::cout << "Bary to Sphere" << std::endl;

    for(std::string &s : subjects)
      c3::bary_to_sphere( params, s );

    return 0;
  }



  //
  // Actual Subj stuff
  //
  c3::Subject_Factory sf;
  std::shared_ptr<c3::Subject> s1 = sf.make_subject(command);
  s1->set(params);
  double* kern_lookup = NULL;
  double** harm_lookup = NULL;

  if(s1->requires_harm_lookup_table()){
    if(verbose > 0)
      std::cout << "Setting Harmonic Lookup" << std::endl;
    harm_lookup = s1->calc_harm_lookup_table(
      c3::string_to_int(params["OPT_VAL_exp_num_harm_samps"]),
      c3::string_to_int(params["OPT_VAL_num_harm"]));
  }

  if(s1->requires_kern_lookup_table()){
    if(verbose > 0)
      std::cout << "Setting Kernel Lookup" << std::endl;

    if (params.find("sigma") == params.end()){
      std::cout << "[c3/main] "
        << "Command " << command << " requires '--sigma' to be set"
        << "Aborting." << std::endl;
      exit(1);
    }
    double sigma = c3::string_to_double(params["sigma"]);
    kern_lookup = s1->calc_kern_lookup_table(sigma,
      c3::string_to_int(params["OPT_VAL_exp_num_kern_samps"]));
  }

  num_subjects = subjects.size();

  //TODO: This is where we'd parallel.
  for(int i = 0; i < num_subjects; ++i){
    std::shared_ptr<c3::Subject> s;
    std::string filename;
    if(i == 0){
      s = s1;
    } else {
      if(s->requires_harm_lookup_table())
        s->set_harm_lookup_table(harm_lookup);
      if(s->requires_kern_lookup_table())
        s->set_kern_lookup_table(kern_lookup);
    }
      

    if(s->requires_trk_endpoints()){
      filename = params["LOAD_xing_path"] + subjects[i] +
          params["LOAD_xing_postfix"];

      if(verbose > 0)
        std::cout << "Loading trk endpoints from " <<
          filename << std::endl;
      s->load_trk_endpoints(filename);
    }
    
    if(s->requires_grid()){
      if(verbose > 0)
        std::cout << "Loading grid" << std::endl;
      s->load_grid(params["LOAD_grid_file"], true, true);
    }

    if(s->requires_kernel()){
      filename = params["LOAD_kernel_path"] + subjects[i] +
        params["LOAD_kernel_postfix"];

      if(verbose > 0)
        std::cout << "Loading kernel" << std::endl;
      s->load_kernel(filename);
    }

    if(s->requires_mask()){
      if(verbose > 0)
        std::cout << "Loading mask" << std::endl;
      s->load_kernel(params["LOAD_mask_file"]);
    }

    int error_code = s->subj_specific_load(subjects[i]);
    if(error_code != 0){
      std::cerr << "Error: s_ops subj_specific_load() error, aborting"
        << std::endl;
      exit(1);
    }

    if(verbose > 0)
      std::cout << "Beginning Run" << std::endl;
    error_code = s->run();
    if(error_code != 0){
      std::cerr << "Error: s_ops run() error, aborting." << std::endl;
      exit(1);
    }

    filename = params["SAVE_" + command + "_prefix"] + subjects[i] +
      params["SAVE_" + command + "_postfix"];
    if(verbose > 0)
      std::cout << "Printing" << std::endl;

    if(!s->no_main_save()){
      s->save_file( filename );
    }

  }

  return 0;
}






