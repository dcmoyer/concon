
#include "trackset.hpp"

namespace trk {

//here, static means that this function does not go outside the TU
static int bytes_to_int (char *b, int offset){
  return *(int*)(&b[offset]) ;
}

static int bytes_to_short (char *b, int offset){
  return *(short*)(&b[offset]) ;
}

Trackset::Trackset(){
  verbose = 0;
  num_tracks = 0;
  max_x = -std::numeric_limits<float>::max();
  max_y = -std::numeric_limits<float>::max();
  max_z = -std::numeric_limits<float>::max();

  min_x = std::numeric_limits<float>::max();
  min_y = std::numeric_limits<float>::max();
  min_z = std::numeric_limits<float>::max();

}

Trackset::Trackset(int _verbose){
  verbose = verbose;
  num_tracks = 0;
  max_x = -std::numeric_limits<float>::max();
  max_y = -std::numeric_limits<float>::max();
  max_z = -std::numeric_limits<float>::max();

  min_x = std::numeric_limits<float>::max();
  min_y = std::numeric_limits<float>::max();
  min_z = std::numeric_limits<float>::max();

}

void Trackset::clear(){

  streams.clear();

  num_tracks = 0;

  max_x = max_y = max_z = 0;
  min_x = min_y = min_z = 0;

}

int Trackset::load_trk(std::string filename){

  std::ifstream trk_file(filename.c_str(), std::ifstream::binary);

  if(!trk_file.good()){
    std::cout << "Bad trk file at " << filename << std::endl;
    return 1;
  }

  trk_file.read(header,1000);

  //TODO:num_tracks verification
  // also when we write the header we should probably change the number of
  // tracks?
  //num_tracks = bytes_to_int_(header,999-12); 

  num_scalars = bytes_to_short(header,36);

  num_props = bytes_to_short(header,236);

  //std::cout << num_tracks << std::endl;

  //return 0;

  int num_points, real_num_tracks = 0;
  //for(int i = 0; i < num_tracks; ++i){
  float coord[3];
  float prop = 0;
  float scalar = 0;

  std::vector<float> scalars;

  trk_file.read((char*)&num_points, sizeof(int));
  while(!trk_file.eof()){

    trk::Streamline s;
    trk::point p;

    scalars.resize( num_scalars * num_points );

    for(int i = 0; i < num_points; ++i){
      trk_file.read((char*)&coord, sizeof(float)*3);

      p.x = coord[0];
      p.y = coord[1];
      p.z = coord[2];
      s.add(p);

      for(int j = 0; j < num_scalars; ++j){
        trk_file.read((char*)&scalar, sizeof(float));
        scalars[i + j*num_points] = scalar;
      }

      if(coord[0] < min_x){
        min_x = coord[0];
      } else if ( coord[0] > max_x ){
        max_x = coord[0];
      }

      if(coord[1] < min_y){
        min_y = coord[1];
      } else if ( coord[1] > max_y ){
        max_y = coord[1];
      }

      if(coord[2] < min_z){
        min_z = coord[2];
      } else if ( coord[2] > max_z ){
        max_z = coord[2];
      }

    }

    for(int i = 0; i < num_props; ++i){
      trk_file.read((char*)&prop, sizeof(float));
      s.add_prop(prop);
    }

    if(num_scalars > 0){
      s.set_all_scalars(scalars);
    }

    streams.push_back(s);
    real_num_tracks += 1;

    trk_file.read((char*)&num_points, sizeof(int));
  }

  trk_file.close();
  //std::cout << real_num_tracks << std::endl;
  num_tracks = real_num_tracks;
  return 0;
}


//
//  save_trk
//
int Trackset::save_trk(std::string filename){

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
  for(int i = 0; i < num_scalars; ++i){
    std::string s_name = "scalar " + std::to_string(i);
    strcpy(&header[38 + 20*i],s_name.c_str());
  }

  //TODO:remove this
  memcpy(&header[236], &num_props, sizeof(short));
  for(int i = 0; i < num_props; ++i){
    std::string s_name = "prop " + std::to_string(i) + "\0";
    strcpy(&header[238 + 20*i],s_name.c_str());
  }

  trk_file.write(header,1000);

  for(int trk_idx = 0; trk_idx < num_tracks; ++trk_idx){

    int num_points;
    float coord[3];
    Streamline& local_stream = streams[trk_idx];

    //bad practice...
    std::vector<float>& local_props = local_stream.props;
    std::vector<float>& local_scalars = local_stream.scalars;

    num_points = local_stream.size();

    trk_file.write(
      reinterpret_cast<char*>(&num_points),
      sizeof(int)
    );

    for(int i = 0; i < num_points; ++i){

      coord[0] = local_stream[i].x;
      coord[1] = local_stream[i].y;
      coord[2] = local_stream[i].z;

      trk_file.write(
        reinterpret_cast<char*>(&coord),
        (sizeof(float)*3)/sizeof(char)
      );

      //
      //  scalars
      //

      if(num_scalars > 0){
        for(int j = 0; j < num_scalars; ++j){
          trk_file.write(
            reinterpret_cast<char*>(&local_scalars[i + j*num_points]),
            (sizeof(float))/sizeof(char)
          );
        }
      }
    }

    //
    //  props
    //

    if(num_props > 0){
      for(int j = 0; j < num_props; ++j){
        trk_file.write(
          reinterpret_cast<char*>(&local_props[j]),
          (sizeof(float)*num_props)/sizeof(char)
        );
      }
    }

  }//end of loop over streams

  trk_file.close();
  return 0;
}


//
//  Prop/Scalar Handlers
//

void Trackset::add_prop(){
  for(Streamline& s : streams){
    s.add_prop(0);
  }
  num_props += 1;
}

void Trackset::add_scalar(){
  for(Streamline& s : streams){
    s.add_scalar( std::vector<float>(s.size(),0) );
  }
  num_scalars += 1;
}

void Trackset::add_prop(const std::vector<float>& values){

  if(values.size() != streams.size()){
    std::cerr << "[trk/trackset:add_prop] " <<
      "values.size() != streams.size(), Aborting." << std::endl;
  }

  for(int i = 0; i < num_tracks; ++i){
    streams[i].add_prop(values[i]);
  }
  num_props += 1;
}

void Trackset::add_scalar(const std::vector<std::vector<float> >& values){

  if(values.size() != streams.size()){
    std::cerr << "[trk/trackset:add_prop] " <<
      "values.size() != streams.size(), Aborting." << std::endl;
  }

  for(int i = 0; i < num_tracks; ++i){
    streams[i].add_scalar(values[i]);
  }
  num_scalars += 1;
}

void Trackset::add_const_scalar(const std::vector<float>& values){

  if(values.size() != streams.size()){
    std::cerr << "[trk/trackset:add_prop] " <<
      "values.size() != streams.size(), Aborting." << std::endl;
  }

  for(int i = 0; i < num_tracks; ++i){
    streams[i].add_scalar( std::vector<float>(streams[i].size(),values[i]) );
  }
  num_scalars += 1;
}

void Trackset::remove_prop(const int& idx){
  for(Streamline &s : streams){
    s.remove_prop(idx);
  }
  num_props -= 1;
}

void Trackset::remove_scalar(const int& idx){
  for(Streamline &s : streams){
    s.remove_scalar(idx);
  }
  num_props -= 1;
}

void Trackset::resample_trajectories(const int& num_points){
  for(Streamline &s : streams)
    s.resample(num_points);
}

//
//  Streamline implementations
//

void Streamline::add_scalar(const std::vector<float>& values){
  if(values.size() != points.size()){
    std::cerr << "[trk/trackset:add_scalar] " <<
      "values size != points size, Aborting." << std::endl;
    exit(1); 
  }
  scalars.insert(scalars.end(), values.begin(), values.end());
}

void Streamline::remove_scalar(const int& idx){
  if(scalars.size() == 0 || (int) (scalars.size() % points.size()) < idx){
    std::cerr << "[trk/trackset:remove_scalar] " <<
      "bad remove idx, Aborting." << std::endl;
    exit(1); 
  }
  scalars.erase(scalars.begin() + points.size() * idx,
    scalars.begin() + points.size() * (idx + 1));
}

//
//  resample
//

void Streamline::resample(const int& num_points){
  if(num_points < 3){
    std::cerr << "[trk/trackset:resample] " <<
      "Resample must have num_poinst > 2. Aborting." << std::endl;
    exit(1);
  }

  int old_num_points = points.size();
  std::vector<point> new_points;
  new_points.resize(num_points);
  new_points[0] = points[0];
  new_points[num_points - 1] = points[num_points - 1];

  //TODO: add support for scalars
  if(scalars.size() > 0){
    std::cerr << "[trk/trackset:resample] " <<
      "Resample + scalars not supported, aborting." << std::endl;
    exit(1);
  }

  int s = num_points - 1;
  int interp_idx;
  float frac;
  for(int i = 1; i < s; ++i){
    interp_idx = (i * (old_num_points - 1)) / (num_points - 1);
    frac = (((float) i * (old_num_points - 1)) / (num_points - 1) )
      - interp_idx;
    new_points[i].x = ((1.0-frac) * points[i].x) + (frac * points[i + 1].x);
    new_points[i].y = ((1.0-frac) * points[i].y) + (frac * points[i + 1].y);
    new_points[i].z = ((1.0-frac) * points[i].z) + (frac * points[i + 1].z);
  }

  points = new_points;

  return;
}


}//end of namespace trk

