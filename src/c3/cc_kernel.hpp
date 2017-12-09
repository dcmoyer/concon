
#ifndef CC_KERNEL_H
#define CC_KERNEL_H

#include <string>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <cassert>
#include <vector>
#include <cmath>

struct cc_kern_triple {
  int row,col;
  double val;
};

struct cc_kern_pair {
  float row,index;
};

class CC_Kernel {
private:
  std::string filename;
  std::vector<cc_kern_triple> rcv;
  std::vector<cc_kern_pair> jump_table;
  int network_size;
  

public:
  CC_Kernel(){;}
  
  int load_file(std::string _filename, bool binary_flag=true);
  int load_binary_file(std::string _filename);
  int load_csv_file(std::string _filename);
  
  double operator()(int row, int col);

  template<typename T>  void to_arma_csc_vectors(std::vector<T> &row,
    std::vector<T> &col, std::vector<double> &value);

  double row_sum(int row);
  int row_num_nonzeros(int row);

  double ordered_array_sum(const int row, const std::vector<int> indices);
  double ordered_array_sum(const int row, const std::vector<int> indices,
  const std::vector<double> node_area);

  double ordered_patch_sum(const std::vector<int> indices_1,
    const std::vector<int> indices_2);
  double ordered_patch_sum(const std::vector<int> indices_1,
    const std::vector<int> indices_2, const std::vector<double> node_area);

  double parc_stats_loglik( const std::vector<int> indices_1,
    const std::vector<int> indices_2, const int mesh_size, const double value);

  double parc_stats_ise( const std::vector<int> indices_1,
    const std::vector<int> indices_2, const int mesh_size);

  double parc_stats_ise( const std::vector<int> indices_1,
    const std::vector<int> indices_2, const int mesh_size, const double average);

  double ordered_array_ise(const int row, const std::vector<int> indices,
    const double value);
  

  void clear(){
    rcv.clear();
    jump_table.clear();
  }
  int size(){
    return network_size;
  }

};

template<typename T>
void CC_Kernel::to_arma_csc_vectors(std::vector<T> &row,
  std::vector<T> &col, std::vector<double> &value){

  row.clear();
  col.clear();
  value.clear();

  int sz = rcv.size();
  for(int i = 0; i < network_size; ++i){
    col.push_back(row.size());
    for(int j = 0; j < sz; ++j){

      if(rcv[j].col < i){
        continue;
      } else if(rcv[j].col == i) {
        row.push_back(rcv[j].row);
        value.push_back(rcv[j].val);
      } else {
        break;
      }

    }
  }
  col.push_back(sz);

}


#endif

