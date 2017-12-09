
#ifndef INTERFACE_HPP
#define INTERFACE_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>

namespace interface {

void default_help();
void default_version();

void load_params(
  int argc, char *argv[],
  std::unordered_map<std::string, std::string>& params,
  std::string& env_file,
  void (*help)() = default_help,
  void (*version)() = default_version,
  int start_index = 1
);

void load_env(
  std::unordered_map<std::string,std::string>& params,
  std::string env_file
);

}// end of namespace interface

#endif



