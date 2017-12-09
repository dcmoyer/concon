
#include "trk.hpp"

namespace trk {

Tractset::Tractset(){
  tracts_x_coords = NULL;
  tracts_x_scalars = NULL;
  tracts_x_props = NULL;

  num_scalars = 0;
  num_props = 0;

  num_points = NULL;
  num_tracts = 0;
  verbose = 0;
}

Tractset::Tractset(int _verbose){
  tracts_x_coords = NULL;
  tracts_x_scalars = NULL;
  tracts_x_props = NULL;

  num_scalars = 0;
  num_props = 0;

  num_points = NULL;
  num_tracts = 0;
  verbose = _verbose;

}

Tractset::~Tractset(){
  clear();
}

int bytes_to_int (char *b, int offset){
  return *(int*)(&b[offset]) ;
}
int bytes_to_short (char *b, int offset){
  return *(short*)(&b[offset]) ;
}

int Tractset::load_trk(std::string filename){

  //filesize stuff
  std::ifstream motd(filename.c_str(),
    std::ifstream::binary|std::ifstream::ate);
  std::ifstream::pos_type size = motd.tellg();
  motd.close();

  std::ifstream trk_file(filename.c_str(), std::ifstream::binary);
  if(!trk_file.good()){
    std::cout << "[trk.cpp:load_trk] Bad trk file at " << filename
      << " aborting." << std::endl;
    return 1;
  }

  //unfortunately we cannot trust this header, because
  //not everyone writes to it e.g. dipy and nibabel
  //char header[1000];
  trk_file.read(header,1000);

  num_scalars = bytes_to_short(header,36);
  num_props = bytes_to_short(header,236);

  //int num_tracks = bytes_to_int(header,1000-4); 
  //std::cout << num_tracks << std::endl;
  //num_tracks = bytes_to_int(header,1000-8); 
  //std::cout << num_tracks << std::endl;
  //num_tracks = bytes_to_int(header,1000-12); 
  //std::cout << num_tracks << std::endl;

  if(verbose > 0){
    std::cout << "[trk] num scalars "<< num_scalars
      << " num properties" << num_props << std::endl;
  }

  //return 0;

  int num_points_local, real_num_tracts = 0;

  //2MB array
  int* num_points_local_array = new int[trk::TRK_LOCAL_ARRAY_SIZE];
  float*** streamline_local_array = new float**[trk::TRK_LOCAL_ARRAY_SIZE];

  float *** tracts_x_scalars_tmp = NULL, ** tracts_x_props_tmp = NULL;

  if(num_scalars > 0)
    tracts_x_scalars_tmp = new float**[trk::TRK_LOCAL_ARRAY_SIZE];
  if(num_props > 0)
    tracts_x_props_tmp = new float*[trk::TRK_LOCAL_ARRAY_SIZE];

  while(trk_file.tellg() < size){
    trk_file.read( reinterpret_cast<char*>(&num_points_local), sizeof(int));

    num_points_local_array[real_num_tracts] = num_points_local;
    streamline_local_array[real_num_tracts] = new float*[num_points_local];

    if(num_scalars > 0)
      tracts_x_scalars_tmp[real_num_tracts] = new float*[num_points_local];
    if(num_props > 0)
      tracts_x_props_tmp[real_num_tracts] = new float[num_props];

    for(int i = 0; i < num_points_local; ++i){
      //streamline[i] = new float[3];
      streamline_local_array[real_num_tracts][i] = new float[3];
      trk_file.read(
        reinterpret_cast<char*>(streamline_local_array[real_num_tracts][i]),
        (sizeof(float)*3)/sizeof(char)
      );

      if(num_scalars > 0){
        tracts_x_scalars_tmp[real_num_tracts][i] = new float[num_scalars];
        trk_file.read(
          reinterpret_cast<char*>(tracts_x_scalars_tmp[real_num_tracts][i]),
          (sizeof(float)*num_scalars)/sizeof(char)
        );
      }

    }

    if(num_props > 0){
      trk_file.read(
        reinterpret_cast<char*>(tracts_x_props_tmp[real_num_tracts]),
        (sizeof(float)*num_props)/sizeof(char)
      );
    }

    real_num_tracts += 1;

    if(real_num_tracts >= trk::TRK_LOCAL_ARRAY_SIZE){
      std::cerr << "...really, 2 million tracts same file?." << std::endl;
      std::cerr << "TODO: expand array here, sorry. Aborting." << std::endl;
      //TODO: expand array here
      exit(1);
    }
  }

  trk_file.close();

  //copy over
  num_tracts = real_num_tracts;
  tracts_x_coords = new float**[num_tracts];
  num_points = new int[num_tracts];

  for(int i = 0; i < num_tracts; ++i){
    //continue;
    tracts_x_coords[i] = streamline_local_array[i];
    num_points[i] = num_points_local_array[i];
  }

  tracts_x_scalars = new float**[num_tracts];
  if(num_scalars > 0){
    for(int i = 0; i < num_tracts; ++i)
      tracts_x_scalars[i] = tracts_x_scalars_tmp[i];
    delete[] tracts_x_scalars_tmp;
  }

  tracts_x_props = new float*[num_tracts];
  if(num_props > 0){
    for(int i = 0; i < num_tracts; ++i)
      tracts_x_props[i] = tracts_x_props_tmp[i];
    delete[] tracts_x_props_tmp;
  }

  delete[] streamline_local_array;
  delete[] num_points_local_array;

  return 0;
}


//
//  save_trk
//
int Tractset::save_trk(std::string filename){

  std::ofstream trk_file(filename.c_str(), std::ifstream::binary);
  if(!trk_file.good()){
    std::cout << "Bad (output) trk file at " << filename << std::endl;
    return 1;
  }

  //unfortunately we cannot trust this header, because
  //not everyone writes to it e.g. dipy and nibabel
  //char header[1000];
  //num_scalars = bytes_to_short(header,36);
  //num_props = bytes_to_short(header,236);

  //TODO:remove this
  memcpy(&header[36], &num_scalars, sizeof(short));
  //short* n_s = (short*)(&header[36]);
  //*n_s = num_scalars;
  for(int i = 0; i < num_scalars; ++i){
    std::string s_name = "scalar " + std::to_string(i);
    strcpy(&header[38 + 20*i],s_name.c_str());
  }

  //TODO:remove this
  memcpy(&header[236], &num_props, sizeof(short));
  //short* n_p = (short*)(&header[236]);
  //*n_p = num_props;
  for(int i = 0; i < num_props; ++i){
    std::string s_name = "prop " + std::to_string(i) + "\0";
    strcpy(&header[238 + 20*i],s_name.c_str());
  }

  trk_file.write(header,1000);

  for(int trk_idx = 0; trk_idx < num_tracts; ++trk_idx){
    trk_file.write(
      reinterpret_cast<char*>(&num_points[trk_idx]),
      sizeof(int)
    );

    int s = num_points[trk_idx];
    for(int i = 0; i < s; ++i){
      trk_file.write(
        reinterpret_cast<char*>(tracts_x_coords[trk_idx][i]),
        (sizeof(float)*3)/sizeof(char)
      );

      if(num_scalars > 0){
        trk_file.write(
          reinterpret_cast<char*>(tracts_x_scalars[trk_idx][i]),
          (sizeof(float)*num_scalars)/sizeof(char)
        );
      }
    }

    if(num_props > 0){
      trk_file.write(
        reinterpret_cast<char*>(tracts_x_props[trk_idx]),
        (sizeof(float)*num_props)/sizeof(char)
      );
    }
  }

  trk_file.close();
  return 0;
}




//
//  clear
//    -deletes current data (and frees memory), but does not dealloc instance
//
void Tractset::clear(){
  
  //TODO: exceptions
  if(tracts_x_coords == NULL && num_points == NULL){
    return;
  } else if (tracts_x_coords == NULL){
    std::cerr << "Bad tracts_x_coords pointer in TKS destructor, Aborting.";
    exit(1);
  } else if (num_points == NULL){
    std::cerr << "Bad num_points pointer in TKS destructor, Aborting.";
    exit(1);
  }  

  for(int i = 0; i < num_tracts; ++i){
    //int s = get_num_points(i);
    int s = num_points[i];
    for(int j = 0; j < s; ++j){
      delete[] tracts_x_coords[i][j];
    }
    delete[] tracts_x_coords[i];
  }

  delete[] tracts_x_coords;
  tracts_x_coords = NULL;

  //
  // dealloc scalars
  //
  if(num_scalars > 0){

    for(int i = 0; i < num_tracts; ++i){
      int s = num_points[i];
      for(int j = 0; j < s; ++j){
        delete[] tracts_x_scalars[i][j];
      }
      delete[] tracts_x_scalars[i];
    }
  }

  delete[] tracts_x_scalars;
  tracts_x_scalars = NULL;
  num_scalars = 0;

  //
  // dealloc props
  //
  if(num_props > 0){
    for(int i = 0; i < num_tracts; ++i){
      delete[] tracts_x_props[i];
    }

  }

  delete[] tracts_x_props;
  tracts_x_props = NULL;
  num_props = 0;

  //
  // dealloc num_points
  //

  delete[] num_points;
  num_points = NULL;

  num_tracts = 0;

}


//
//  scalars and props stuff
//

void Tractset::set_scalar(
  const int& trk_idx, const int& pt_idx, const int& s_idx,
  const float& value){
  if(trk_idx > num_tracts || trk_idx < 0){
    std::cerr << "[struct/trk.cpp:set_scalar] " <<
      "trk_idx bad. Aborting." << std::endl;
    exit(1);
  }

  if(pt_idx > num_points[trk_idx] || pt_idx < 0){
    std::cerr << "[struct/trk.cpp:set_scalar] " <<
      "pt_idx bad. Aborting." << std::endl;
    exit(1);
  }

  if(s_idx > num_scalars || s_idx < 0){
    std::cerr << "[struct/trk.cpp:set_scalar] " <<
      "s_idx bad. Aborting." << std::endl;
    exit(1);
  }

  tracts_x_scalars[trk_idx][pt_idx][s_idx] = value;
}

void Tractset::set_prop(const int& trk_idx, const int& p_idx,
  const float& value){
  if(trk_idx > num_tracts || trk_idx < 0){
    std::cerr << "[struct/trk.cpp:set_prop] " <<
      "trk_idx bad. Aborting." << std::endl;
    exit(1);
  }

  if(p_idx > num_props || p_idx < 0){
    std::cerr << "[struct/trk.cpp:set_prop] " <<
      "p_idx bad. Aborting." << std::endl;
    exit(1);
  }

  tracts_x_props[trk_idx][p_idx] = value; 
}

void Tractset::add_scalar(){
  float ** temp_scalar;
  for(int i = 0; i < num_tracts; ++i){
    temp_scalar = new float*[num_points[i]];

    for(int j = 0; j < num_points[i]; ++j){
      temp_scalar[j] = new float[num_scalars + 1];

      for(int k = 0; k < num_scalars; ++k)
        temp_scalar[j][k] = tracts_x_scalars[i][j][k];
      temp_scalar[j][num_scalars] = 0;

      if(num_scalars > 0)
        delete[] tracts_x_scalars[i][j];
    }

    if(num_scalars > 0)
      delete[] tracts_x_scalars[i];
    tracts_x_scalars[i] = temp_scalar;
      
  }
  num_scalars += 1;


  return;
}

void Tractset::add_prop(){
  float * temp_prop;
  for(int i = 0; i < num_tracts; ++i){
    temp_prop = new float[num_props + 1];
    for(int j = 0; j < num_props; ++j)
      temp_prop[j] = tracts_x_props[i][j];
    temp_prop[num_props] = 0;

    if(num_props > 0)
      delete[] tracts_x_props[i];
    tracts_x_props[i] = temp_prop;
  }

  num_props += 1;
}


//void Tractset::add_scalar(const std::vector<float>& scalar);

//this is a dumb hack because trackvis doesn't support its own format
void Tractset::add_trk_const_scalar(std::vector<float> scalar){
  if((int) scalar.size() != num_tracts){
    std::cerr << "[struct/trk.cpp:add_prop] " <<
      "prop vec size not equal to num_tracts. Aborting." << std::endl;
    exit(1);
  }

  add_scalar();
  int s;
  for(int i = 0; i < num_tracts; ++i){
    s = num_points[i];
    for(int j = 0; j < s; ++j)
      tracts_x_scalars[i][j][num_scalars - 1] = scalar[i];
  }

}


void Tractset::add_prop(const std::vector<float>& prop){
  if((int) prop.size() != num_tracts){
    std::cerr << "[struct/trk.cpp:add_prop] " <<
      "prop vec size not equal to num_tracts. Aborting." << std::endl;
    exit(1);
  }

  add_prop();
  for(int i = 0; i < num_tracts; ++i){
    tracts_x_props[i][num_props - 1] = prop[i];
  }
}

void Tractset::remove_scalar(const int& idx){
  if(idx >= num_scalars || idx < 0){
    std::cerr << "[struckt/trk.cpp:remove_scalar] " <<
      "bad idx. Aborting." << std::endl;

  }
  float ** temp_scalar;
  int s;
  for(int i = 0; i < num_tracts; ++i){

    s = num_points[i];
    temp_scalar = new float*[s];
    for(int j = 0; j < s; ++j){

      if(num_scalars > 1){
        temp_scalar[j] = new float[num_scalars - 1];

        for(int k = 0; k < idx; ++k)
          temp_scalar[j][k] = tracts_x_scalars[i][j][k];
        for(int k = idx + 1; k < num_scalars; ++k)
          temp_scalar[j][k - 1] = tracts_x_scalars[i][j][k];
      }
      delete[] tracts_x_scalars[i][j];

    }
    delete[] tracts_x_scalars[i];
    if(num_scalars > 1){
      tracts_x_scalars[i] = temp_scalar;
    } else {
      delete[] temp_scalar;
    }
  }

  num_scalars -= 1;
  return;
}

void Tractset::remove_prop(const int& idx){
  if(idx >= num_props || idx < 0){
    std::cerr << "[struct/trk.cpp:remove_prop] " <<
      "bad idx. Aborting." << std::endl;
    exit(1);
  }

  float * temp_prop;
  for(int i = 0; i < num_tracts; ++i){

    if(num_props > 1){
      temp_prop = new float[num_props - 1];
      for(int j = 0; j < idx; ++j)
        temp_prop[j] = tracts_x_props[i][j];
      for(int j = idx + 1; j < num_props; ++j)
        temp_prop[j - 1] = tracts_x_props[i][j];
      delete[] tracts_x_props[i];
      tracts_x_props[i] = temp_prop;
    } else {
      delete[] tracts_x_props[i];
    } 

  }

  num_props -= 1;
  return;
}




//
//  Resamples tracts to a regular length
//
void Tractset::resample_trajectories(int N){

  if(tracts_x_coords == NULL || num_points == NULL){
    std::cerr << "Bad call to Trackset::resample, Aborting." << std::endl;
    exit(1);
  }

  if(N < 3){
    std::cerr << "Resample must have N>2. Aborting." << std::endl;
    exit(1);
  }

  //TODO: this.
  if(num_scalars > 0){
    std::cerr << "[struct/trk.cpp] " <<
      "Resample + scalars not supported, aborting." << std::endl;
    exit(1);
  }

  float** tmp_array;
  float** curr_trk; //annoying, yes. Sorry.
  float frac;
  int old_num_points;
  int interp_idx;
  int s;

  if(verbose > 0){
    std::cout << "[struct/trk.cpp] resample trajectories starting"
      << std::endl;
  }

  for(int trk_idx = 0; trk_idx < num_tracts; ++trk_idx){
    if(verbose > 0){
      std::cout << "\r" << trk_idx << std::flush;
    }

    tmp_array = new float*[N];
    old_num_points = num_points[trk_idx];
    curr_trk = tracts_x_coords[trk_idx];

    //check if ill-formed
    if(old_num_points < 0){
      std::cerr << "[struct/trk.cpp] trk_idx " << trk_idx
        << " has no points. Aborting."
        << std::endl;
      exit(1);
    }

    //check if single point
    s = N - 1;
    if(old_num_points == 1){

      //if(verbose > 0){
      //  std::cout << "[trk] single point?" << std::endl;
      //}

      for(int i = 0; i < s; ++i){
        tmp_array[i] = new float[3];
        tmp_array[i][0] = curr_trk[0][0];
        tmp_array[i][1] = curr_trk[0][1];
        tmp_array[i][2] = curr_trk[0][2];
      }
      s = 0; //skips next loop
    }

    //we have N - 1 points instead of old_num_points - 1
    //the - 1 comes from the 0th point being the same
    //However, the last point will also be the same.
    //This would usually work out correctly, except we also might read off
    //the end of the array. So to avoid that, we loop to N - 2 (inclusive)
    //and then tag on the last point at the end.
    //at the (N-1)th step we'd have
    //  frac = 0
    //  interp_idx = old_num_points - 1
    //but then
    //  interp_idx + 1 > array last idx
    //and this causes a seg fault
    for(int i = 0; i < s; ++i){
      interp_idx = (i * (old_num_points - 1)) / (N - 1);

      frac = (((float) i * (old_num_points - 1)) / (N - 1) ) - interp_idx;

      tmp_array[i] = new float[3];
      tmp_array[i][0] = ((1.0 - frac) * curr_trk[interp_idx][0])
        + (frac * curr_trk[interp_idx + 1][0]);
      tmp_array[i][1] = ((1.0 - frac) * curr_trk[interp_idx][1])
        + (frac * curr_trk[interp_idx + 1][1]);
      tmp_array[i][2] = ((1.0 - frac) * curr_trk[interp_idx][2])
        + (frac * curr_trk[interp_idx + 1][2]);

    }
    tmp_array[N - 1] = new float[3];
    tmp_array[N - 1][0] = curr_trk[old_num_points - 1][0];
    tmp_array[N - 1][1] = curr_trk[old_num_points - 1][1];
    tmp_array[N - 1][2] = curr_trk[old_num_points - 1][2];

    //delete old coords
    for(int i = 0; i < old_num_points; ++i){
      delete[] tracts_x_coords[trk_idx][i];
    }
    //delete old array of coords
    delete[] tracts_x_coords[trk_idx];

    //set new array
    tracts_x_coords[trk_idx] = tmp_array;
    num_points[trk_idx] = N;

    tmp_array = NULL;

  }

  if(verbose > 0){
    std::cout << std::endl << 
      "[struct/trk.cpp] resample trajectories exiting" << std::endl;
  }
  return;
}


//
//
// Inefficient! Avoid!
//
void Tractset::delete_tract(int idx){
  if(idx > num_tracts || idx < 0){
    std::cerr << "[struct/trk.cpp:delete_tract] " << 
      "Bad idx. Aborting." << std::endl;
    exit(1);
  }

  //delete coordinates
  for(int i = 0; i < num_points[i]; ++i){
    delete[] tracts_x_coords[idx][i];
    if(num_scalars > 0)
      delete[] tracts_x_scalars[idx][i];
  }

  //delete row of coordinate
  delete[] tracts_x_coords[idx];
  if(num_scalars > 0)
    delete[] tracts_x_scalars[idx];

  if(num_props > 0)
    delete[] tracts_x_props[idx];

  //move everyone forward one!
  for(int i = idx + 1; i < num_tracts; ++i){
    tracts_x_coords[i - 1] = tracts_x_coords[i];
    num_points[i - 1] = num_points[i];

    if(num_scalars > 0)
      tracts_x_scalars[i - 1] = tracts_x_scalars[i];

    if(num_props > 0)
      tracts_x_props[i - 1] = tracts_x_props[i];
  }

  num_tracts -= 1;
  tracts_x_coords[num_tracts] = NULL;
  tracts_x_scalars[num_tracts] = NULL;
  tracts_x_props[num_tracts] = NULL;
  return;
}


void Tractset::filter_by_const_scalar(int scalar_idx,
  std::function<bool (float)> condition){
  //AVOID THE DELETE METHOD
  //
  //

  if(scalar_idx >= num_scalars || scalar_idx < 0){
    std::cerr << "[struct/trk.cpp:filter_by_const_scalar] "
      << "Bad scalar_idx, aborting." << std::endl;
    exit(1);
  }

  int num_to_keep = 0;

  for(int i = 0; i < num_tracts; ++i){
    if( condition(tracts_x_scalars[i][0][scalar_idx]) )
      num_to_keep += 1;
  }

  float *** new_txc = new float**[num_to_keep];
  float *** new_txs = new float**[num_to_keep];
  float ** new_txp = new float*[num_to_keep];
  int * new_num_points = new int[num_to_keep];

  int new_idx = 0;
  for(int i = 0; i < num_tracts; ++i){
    if( condition(tracts_x_scalars[i][0][scalar_idx]) ){
      //copy over
      new_txc[new_idx] = tracts_x_coords[i];
      new_num_points[new_idx] = num_points[i];

      if(num_scalars > 0)
        new_txs[new_idx] = tracts_x_scalars[i];

      if(num_props > 0)
        new_txp[new_idx] = tracts_x_props[i];

      new_idx += 1;
    } else {

      //delete
      for(int j = 0; j < num_points[i]; ++j){
        delete[] tracts_x_coords[i][j];

        if(num_scalars > 0)
          delete[] tracts_x_scalars[i][j];
      }
      delete[] tracts_x_coords[i];

      if(num_scalars > 0)
        delete[] tracts_x_scalars[i];

      if(num_props > 0)
        delete[] tracts_x_props[i];

    }
  }

  delete[] tracts_x_coords;
  delete[] num_points;
  delete[] tracts_x_scalars;
  delete[] tracts_x_props;

  tracts_x_coords = new_txc;
  tracts_x_scalars = new_txs;
  tracts_x_props = new_txp;
  num_points = new_num_points;
  num_tracts = num_to_keep;

  return;
}


void Tractset::filter_by_prop(int prop_idx,
  std::function<bool (float)> condition){
  //AVOID THE DELETE METHOD
  //
  //

  if(prop_idx > num_props || prop_idx < 0){
    std::cerr << "[struct/trk.cpp:filter_by_prop] "
      << "Bad prop_idx, aborting." << std::endl;
    exit(1);
  }

  int num_to_keep = 0;

  for(int i = 0; i < num_tracts; ++i){
    if( condition(tracts_x_props[i][prop_idx]) )
      num_to_keep += 1;
  }

  float *** new_txc = new float**[num_to_keep];
  float *** new_txs = new float**[num_to_keep];
  float ** new_txp = new float*[num_to_keep];
  int * new_num_points = new int[num_to_keep];

  int new_idx = 0;
  for(int i = 0; i < num_tracts; ++i){
    if( condition(tracts_x_props[i][prop_idx]) ){
      //copy over
      new_txc[new_idx] = tracts_x_coords[i];
      new_num_points[new_idx] = num_points[i];

      if(num_scalars > 0)
        new_txs[new_idx] = tracts_x_scalars[i];

      if(num_props > 0)
        new_txp[new_idx] = tracts_x_props[i];

      new_idx += 1;
    } else {

      //delete
      for(int j = 0; j < num_points[i]; ++j){
        delete[] tracts_x_coords[i][j];

        if(num_scalars > 0)
          delete[] tracts_x_scalars[i][j];
      }
      delete[] tracts_x_coords[i];

      if(num_scalars > 0)
        delete[] tracts_x_scalars[i];

      if(num_props > 0)
        delete[] tracts_x_props[i];

    }
  }

  delete[] tracts_x_coords;
  delete[] num_points;
  delete[] tracts_x_scalars;
  delete[] tracts_x_props;

  tracts_x_coords = new_txc;
  tracts_x_scalars = new_txs;
  tracts_x_props = new_txp;
  num_points = new_num_points;
  num_tracts = num_to_keep;

  return;
}





//
//  resample trajectories
//

void Tractset::resample_tracts(int N, bool replacement, bool in_order,
  int seed){

  std::srand(seed);

  if(tracts_x_coords == NULL || num_points == NULL){
    std::cerr << "[struct/trk.cpp:resample_tracts]" <<
      "Bad call to Trackset::resample, Aborting." << std::endl;
    exit(1);
  }

  if(N < 0){
    std::cerr << "[struct/trk.cpp:resample_tracts]" <<
     "Resample must have N>0. Aborting." << std::endl;
    exit(1);
  }

  if(verbose > 0){
    std::cout << "[struct/trk.cpp] resample tracts starting" << std::endl;
  }

  if(replacement && in_order){
    std::cerr << "[struct/trk.cpp:resample_tracts] replacement and in_order"
      << " specified. Aborting." << std::endl;
    exit(1);
  }

  if(!replacement && N > num_tracts){
    std::cerr << "[struct/trk.cpp:resample_tracts] no replacement and "
      << "N > num_tracts. Aborting." << std::endl;
    exit(1);
  }

  //TODO: this.
  if(num_props > 0){
    std::cerr << "[struct/trk.cpp] " <<
      "Resample + props not supported, aborting." << std::endl;
    exit(1);
  }

  std::vector<int> possible_indices;
  possible_indices.reserve(num_tracts);
  for(int i = 0; i < num_tracts; ++i){
    possible_indices.push_back(i);
  }

  std::vector<int> resample_indices;
  resample_indices.reserve(N);

  int idx;
  for(int i = 0; i < N; ++i){

    if(replacement){

      idx = std::rand() % num_tracts;
      resample_indices.push_back(idx);
      possible_indices[idx] = -1;

    } else {

      do {
        idx = std::rand() % num_tracts;
      } while(possible_indices[idx] == -1);
      resample_indices.push_back(idx);
      possible_indices[idx] = -1;

    }
  }

  if(in_order){
    std::sort(resample_indices.begin(), resample_indices.end());
  }

  float *** new_trk_array = new float**[N];
  int * new_num_points = new int[N];
  for(int trk_idx = 0; trk_idx < N; ++trk_idx){
    new_trk_array[trk_idx] = tracts_x_coords[resample_indices[trk_idx]];
    new_num_points[trk_idx] = num_points[resample_indices[trk_idx]];
  }

  //deallocation
  for(int i = 0; i < num_tracts; ++i){

    if(possible_indices[i] > -1){
      for(int s = 0; s < num_points[i]; ++s){
        delete[] tracts_x_coords[i][s];
      }
      delete[] tracts_x_coords[i];

    }
  }
  delete[] num_points;
  delete[] tracts_x_coords;

  num_points = new_num_points;
  tracts_x_coords = new_trk_array;
  num_tracts = N;

  if(verbose > 0){
    std::cout << std::endl << 
      "[struct/trk.cpp] resample tracts exiting" << std::endl;
  }
  return;
}





}//end of namespace trk

