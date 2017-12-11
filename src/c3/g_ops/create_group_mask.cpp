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



#include "c3/g_ops/create_group_mask.hpp"



int c3::create_group_mask( std::vector<std::string> filenames,
  std::string output, bool verbose){

  //check if file exists
  for( std::string file : filenames ){
    std::ifstream check(file.c_str()); 
    if( !check.good()){
      std::cerr << "Bad Filename List at file: " << file << "Aborting"
        << std::endl;
      return 1;
    }
    check.close();
  }

  int network_size = -1;
  int row, col;
  double value;
  bool ** mask = NULL;
  //load each file in order
  for( std::string file : filenames ){

    if(verbose)
      std::cout << file << std::endl;

    std::ifstream reader(file.c_str(), std::ios::binary);

    //check network size
    if(network_size == -1) {
      reader.read(reinterpret_cast<char*>(&network_size),sizeof(int));
    } else {
      int temp_network_size;
      reader.read(reinterpret_cast<char*>(&temp_network_size),sizeof(int));
      if(temp_network_size != network_size){
        std::cerr << "Bad network size at file: " << file << "Aborting"
          << std::endl;
        return 1;
      }
    }

    //set up mask
    if(mask == NULL){
      mask = new bool*[network_size];
      for(int i = 0; i < network_size; ++i){
        mask[i] = new bool[network_size];
        for(int j = 0; j < network_size; ++j)
          mask[i][j] = false;
      }
    }

    //fill mask
    while( !reader.eof() ){

      reader.read( reinterpret_cast<char*>(&row),sizeof(int));
      reader.read( reinterpret_cast<char*>(&col),sizeof(int));
      reader.read( reinterpret_cast<char*>(&value),sizeof(double));
      mask[row][col] = true;
      //if(verbose){
      //  std::cout << "\rR" << row 
      //            << "\tC" << col
      //            << "\tV"<< value << std::flush;
      //}

    }
    if(verbose)
      std::cout << std::endl;

    reader.close();
  }

  std::ofstream writer(output.c_str());
  writer << network_size << std::endl;
  for(int i = 0; i < network_size; ++i)
    for(int j = i; j < network_size; ++j)
      if(mask[i][j])
        writer << i << "\t" << j << std::endl;
  writer.close();

  return 0;
}


