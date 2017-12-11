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


#include "c3/g_ops/merge_marginals.hpp"


namespace c3 {

int merge_marginals( std::vector<std::string> subj_list,
  std::unordered_map<std::string,std::string> params){

  std::string out_prefix, out_postfix;
  bool verbose = false;
  std::vector<int> target_regions;

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

    if(p.first == "out_prefix"){
      out_prefix = p.second;
    }

    if(p.first == "out_postfix"){
      out_postfix = p.second;
    }
  }

  std::string region_string;
  for( int& region : target_regions ){

    if(region == -2){
      region_string = "full";
    } else {
      region_string = std::to_string(region);
    }

    std::ofstream agg_output( out_prefix + "agg_" + region_string +
      out_postfix);
    bool first = true;
    for( std::string &subj : subj_list ){

      if(!first){
        agg_output << std::endl;
      } else {
        first = false;
      }

      std::ifstream in_subj( out_prefix + subj + "_" + region_string
        + out_postfix);

      if(!in_subj.good()){
        in_subj.close();
        first = true;
        continue;
      }

      agg_output << subj;
      double temp;
      while(in_subj >> temp){
        agg_output << "\t" << temp;
      }
      in_subj.close();

    }
    agg_output.close();
  }

  return 0;
}





}// end of namespace c3

