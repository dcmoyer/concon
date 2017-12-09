
#include "interface.hpp"

namespace interface {

void default_help(){
  std::cout << "Default Interface Help, triggered by -h arg pass."
    << std::endl;
  return;
}

void default_version(){
  std::cout << "Default Interface Version, triggered by --version arg pass."
    << std::endl;
  return;
}


void load_params(
  int argc, char *argv[],
  std::unordered_map<std::string, std::string>& params,
  std::string& env_file,
  void (*help)(),
  void (*version)(),
  int start_index 
){

  //put args into params
  std::string arg,value;
  int verbose = 0;

  //TODO: single arg configuration?
  for(int i = start_index; i < argc; i += 2){

    arg = std::string(argv[i]);

    //help function
    if(arg == "-h" || arg == "--help"){
      help();
      exit(0);
    }

    //help function
    if(arg == "--version"){
      version();
      exit(0);
    }

    //TODO: single flag args?
    if(argc <= i+1){
      std::cout << "[interface] incorrect number of args" << std::endl;
      exit(1);
    }
    value = std::string(argv[i+1]);    

    //TODO: generalize?
    if(arg.length() < 3){
      std::cout << "[interface] arg 2char too short" << std::endl;
      exit(1);
    }

    //TODO: generalize???
    if(arg.substr(0,2) != "--"){
      std::cout << "[interface] arg not -- at the start" << std::endl;
      exit(1);
    }
    arg.erase(0,2);

    //replace all -'s with _'s
    int s = arg.length();
    for(int idx = 0; idx < s; ++idx){
      if(arg[idx] == '-'){
        arg[idx] = '_';
      }
    }

    if(arg == "verbose"){
      verbose = 1;
    }

    if(verbose > 0){
      std::cout << "[interface] INPUT: " << arg << " " << value << std::endl;
    }

    if(arg == "env" || arg == "env_file"){
      env_file = value;
      continue;
    }

    params[arg] = std::string(value);
  }
}

void load_env(
  std::unordered_map<std::string,std::string>& params,
  std::string env_file
){

  std::string input_key, input_value;
  int verbose;
  std::ifstream env_input(env_file.c_str());

  while( env_input >> input_key >> input_value ){
    if(input_key == "verbose"){
      verbose = 1;
    }

    if(verbose > 0){
      std::cout << "[interface] INPUT (env-file): " << input_key <<
        " " << input_value << std::endl;
    }

    if(params.find(input_key) == params.end()){
      params[input_key] = input_value;
    } else {

      if(verbose > 0)
        std::cout << "[interface]: " << input_key <<
          " overriden by cmd line." << std::endl;

    }

  }
  env_input.close();

}


}//end of namespace interface


