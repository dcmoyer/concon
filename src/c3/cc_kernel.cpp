#include "c3/cc_kernel.hpp"


int CC_Kernel::load_file(std::string _filename, bool binary_flag) {
  
  if(binary_flag){
    return load_binary_file( _filename);
  } else {
    return load_csv_file( _filename);
  } 
}

int CC_Kernel::load_binary_file(std::string _filename){ 
  int counter,previous_row;
  FILE *ptr_myfile = NULL;
  struct cc_kern_triple row_col_val_triple;
  struct cc_kern_pair row_index_pair;
  filename = _filename;

  ptr_myfile = fopen(filename.c_str(), "rb");
  if (!ptr_myfile) {
    //printf("Unable to open file!\n");
    std::cout << "Unable to open file " << filename << std::endl
      << "Aborting." << std::endl;
    assert(!!ptr_myfile);
    return 1;
  }

  previous_row = -1;
  counter = -1;
  //network_size = 0;
  fread( &network_size, sizeof(int), 1, ptr_myfile);
  int row, col;
  double val;
  while( true ) {
    
    //std::cout<< counter << std::endl;
    //fread( &row_col_val_triple, sizeof(struct cc_kern_triple), 1, ptr_myfile);

    fread( &row, sizeof(int), 1, ptr_myfile);
    fread( &col, sizeof(int), 1, ptr_myfile);
    fread( &val, sizeof(double), 1, ptr_myfile);
    row_col_val_triple.row = row;
    row_col_val_triple.col = col;
    row_col_val_triple.val = val;

    //std::cout << row << "\t";
    //std::cout << col << "\t";
    //std::cout << val << std::endl;
 
    //std::cout << "WARNING, ASSUMING 1 INDEXING USING BORIS' SPARSE MATRICES";
    //row_col_val_triple.row -= 1;
    //row_col_val_triple.col -= 1;

    if(!feof(ptr_myfile)){
      rcv.push_back(row_col_val_triple);

      //if(row_col_val_triple.col > network_size)
      //  network_size = row_col_val_triple.col;
      //if(row_col_val_triple.row > network_size)
      //  network_size = row_col_val_triple.row;

      counter += 1;
      if(row_col_val_triple.row > previous_row){
        previous_row = row_col_val_triple.row;
        row_index_pair.row = previous_row;
        row_index_pair.index = counter;
        jump_table.push_back(row_index_pair);
      }
    } else {
      break;
    }
  }
  
  fclose(ptr_myfile);
  return 0;

}

int CC_Kernel::load_csv_file( std::string _filename){
  std::cout << "loading: " << _filename << std::endl;
  std::ifstream s( _filename.c_str(), std::ifstream::in);
  double temp;
  char c;
 
  int i = 0, j = 0, max_rows = 1, counter = 0, previous_row = -1;
  struct cc_kern_triple row_col_val_triple;
  struct cc_kern_pair row_index_pair;
  s >> network_size;
  while( i < max_rows ){
   
    s >> temp;
    if(temp != 0){
      row_col_val_triple.row = i;
      row_col_val_triple.col = j;
      row_col_val_triple.val = temp;
      rcv.push_back(row_col_val_triple);

      if(i > previous_row){
        previous_row = i;
        row_index_pair.row = previous_row;
        row_index_pair.index = counter;
        jump_table.push_back(row_index_pair);
      }

      counter += 1;
    }
    
    s.get(c);
    //std::cout << c;
    if( c == ',' || c == '\t'){
      j += 1;
    } else {
      max_rows = j + 1;
      j = 0;
      i += 1;
    }
  }
  s.close();
  //network_size = max_rows;
  return 0;
}

double CC_Kernel::operator() (int row, int col){
  int jump = -1;
  int stop = -1;
  int sz = jump_table.size();
  for(int i = 0; i < sz; ++i){
    if(jump_table[i].row == row){
      jump = jump_table[i].index;
      if(i < sz - 1)
        stop = jump_table[i+1].index;
      else
        stop = rcv.size();
      break;
    }
  }

  //std::cout << jump << std::endl;
      
  if(jump == -1){
    return 0.0;
  }

  for(int j = jump; j < stop ; ++j){
    if(rcv[j].col == col){
      return rcv[j].val;
    } else if (rcv[j].col > col){
      break;
    }
  }
  return 0.0;
}

double CC_Kernel::row_sum(int row){
  int jump = -1;
  int stop = -1;
  int sz = jump_table.size();
  for(int i = 0; i < sz; ++i){
    if(jump_table[i].row == row){
      jump = jump_table[i].index;
      if(i < sz - 1)
        stop = jump_table[i+1].index;
      else
        stop = rcv.size();
      break;
    }
  }

  //std::cout << jump << std::endl;
      
  if(jump == -1){
    return 0.0;
  }

  double sum = 0;
  for(int j = jump; j < stop ; ++j){
    sum += rcv[j].val;
  }
  return sum;
}

int CC_Kernel::row_num_nonzeros(int row){
  int jump = -1;
  int stop = -1;
  int sz = jump_table.size();
  for(int i = 0; i < sz; ++i){
    if(jump_table[i].row == row){
      jump = jump_table[i].index;
      if(i < sz - 1)
        stop = jump_table[i+1].index;
      else
        stop = rcv.size();
      break;
    }
  }

  //std::cout << jump << std::endl;
      
  if(jump == -1){
    return 0.0;
  }

  int sum = 0;
  for(int j = jump; j < stop ; ++j){
    if(rcv[j].val > 0)
      sum += 1;
  }
  return sum;
}

//
//  Sum of given indices for a specific row.
//  Assumes indices are in order.
//
double CC_Kernel::ordered_array_sum(const int row, const std::vector<int> indices){
  int input_size = indices.size();

  if(input_size == 0){
    return 0.0;
  }

  int jump = -1;
  int stop = -1;
  int sz = jump_table.size();
  for(int i = 0; i < sz; ++i){
    if(jump_table[i].row == row){
      jump = jump_table[i].index;
      if(i < sz - 1)
        stop = jump_table[i+1].index;
      else
        stop = rcv.size();
      break;
    }
  }

  //std::cout << jump << std::endl;
      
  if(jump == -1){
    return 0.0;
  }

  double sum = 0;
  int current_idx = 0;
  for(int j = jump; j < stop ; ++j){
    while(indices[current_idx] < rcv[j].col && current_idx < input_size)
      current_idx++;
    if(current_idx >= input_size)
      break;
    if(indices[current_idx] == rcv[j].col){
      sum += rcv[j].val;
      current_idx++;
    }
  }
  return sum;
}

double CC_Kernel::ordered_array_sum(const int row, const std::vector<int> indices,
  const std::vector<double> node_area){

  int input_size = indices.size();

  if(input_size == 0){
    return 0.0;
  }

  assert((unsigned int) network_size == node_area.size());

  int jump = -1;
  int stop = -1;
  int sz = jump_table.size();
  for(int i = 0; i < sz; ++i){
    if(jump_table[i].row == row){
      jump = jump_table[i].index;
      if(i < sz - 1)
        stop = jump_table[i+1].index;
      else
        stop = rcv.size();
      break;
    }
  }

  //std::cout << jump << std::endl;
      
  if(jump == -1){
    return 0.0;
  }

  double sum = 0;
  int current_idx = 0;
  for(int j = jump; j < stop ; ++j){
    while(indices[current_idx] < rcv[j].col && current_idx < input_size)
      current_idx++;
    if(current_idx >= input_size)
      break;
    if(indices[current_idx] == rcv[j].col){
      sum += rcv[j].val * node_area[rcv[j].col];
      current_idx++;
    }
  }
  return sum;
}

double CC_Kernel::ordered_patch_sum(const std::vector<int> indices_1,
    const std::vector<int> indices_2){

  int size_1 = indices_1.size();
  double output = 0;

  for(int i = 0; i < size_1; ++i){
    output += ordered_array_sum(indices_1[i],indices_2);
  }
  return output;
}

double CC_Kernel::ordered_patch_sum(const std::vector<int> indices_1,
    const std::vector<int> indices_2, const std::vector<double> node_area){

  assert((unsigned int) network_size == node_area.size());
  int size_1 = indices_1.size();
  double output = 0;

  for(int i = 0; i < size_1; ++i){
    output += ordered_array_sum(indices_1[i],indices_2, node_area) *
    node_area[indices_1[i]];
  }
  return output;
}

double CC_Kernel::parc_stats_loglik(const std::vector<int> indices_1,
    const std::vector<int> indices_2, const int mesh_size, const double value){

  double loglik = 0;
  
  int s1 = indices_1.size();
  int s2 = indices_2.size();

  double average = (1.0 / s1) * (1.0 / s2) *
    ordered_patch_sum(indices_1,indices_2);

  //Poisson distribution
  if(average == 0){
    return -20.0;
  }
  loglik += value * log(average) - average;
  for(int i = 1; i < value; ++i)
    loglik += -log(i);

  return fmax(loglik,-20.0);
}

double CC_Kernel::parc_stats_ise( const std::vector<int> indices_1,
    const std::vector<int> indices_2, const int mesh_size){

  int s1 = indices_1.size();
  int s2 = indices_2.size();

  double average = (1.0 / s1) * (1.0 / s2) *
    ordered_patch_sum(indices_1,indices_2);

  double ISE = 0;
  for(int i = 0; i < s1; ++i){
    ISE += (1.0/mesh_size) * (1.0/mesh_size) *
      ordered_array_ise(indices_1[i], indices_2, average);
    //for(int j = 0; j < s2; ++j){
    //  value = (*this)(i,j);
    //  ISE += (1.0/mesh_size) * (1.0/mesh_size) *
    //    (average - value) * (average - value);
    //}
  }

  return ISE;
}

//for computing test-retest ise
double CC_Kernel::parc_stats_ise( const std::vector<int> indices_1,
  const std::vector<int> indices_2, const int mesh_size, const double average){

  int s1 = indices_1.size();

  double ISE = 0;
  for(int i = 0; i < s1; ++i){
    ISE += (1.0/mesh_size) * (1.0/mesh_size) *
      ordered_array_ise(indices_1[i], indices_2, average);
    //for(int j = 0; j < s2; ++j){
    //  value = (*this)(i,j);
    //  ISE += (1.0/mesh_size) * (1.0/mesh_size) *
    //    (average - value) * (average - value);
    //}
  }

  return ISE;
}




double CC_Kernel::ordered_array_ise(const int row, const std::vector<int> indices,
  const double value){

  int input_size = indices.size();

  if(input_size == 0){
    return 0.0;
  }


  int jump = -1;
  int stop = -1;
  int sz = jump_table.size();
  for(int i = 0; i < sz; ++i){
    if(jump_table[i].row == row){
      jump = jump_table[i].index;
      if(i < sz - 1)
        stop = jump_table[i+1].index;
      else
        stop = rcv.size();
      break;
    }
  }

  //std::cout << jump << std::endl;
      
  if(jump == -1){
    return 0.0;
  }

  double sum = 0;
  int current_idx = 0;
  for(int j = jump; j < stop ; ++j){
    while(indices[current_idx] < rcv[j].col && current_idx < input_size)
      current_idx++;
    if(current_idx >= input_size)
      break;
    if(indices[current_idx] == rcv[j].col){
      sum += (rcv[j].val - value) * (rcv[j].val - value);
      current_idx++;
    }
  }
  return sum;
}
