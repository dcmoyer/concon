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



#ifndef TRACKSET_HPP
#define TRACKSET_HPP

#include <string>
#include <string.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
//#include "MeshLib/Point.h"
//TODO: for fast streamline search?
//#include "c3/oct_tree.hpp"

namespace trk {

struct point {
  float x;
  float y;
  float z;
};

//class TrackIterator;
class Streamline;
typedef std::vector<Streamline>::iterator TrackIterator;
//class StreamlineIterator;
typedef std::vector<point>::iterator StreamlineIterator;

class Trackset {
public:
  Trackset();
  Trackset(int _verbose);
  //~Trackset();

  //add streamline?
  void add_streamline(Streamline& s){
    //TODO:update max/mins
    streams.push_back(s);
    num_tracks += 1;
  }

  void delete_streamline(const int& idx){
    streams.erase(streams.begin() + idx);
  }
  //delete streamline?

  int load_ucf(std::string filename);
  int load_trk(std::string filename);

  int save_trk(std::string filename);

  TrackIterator get_iterator(){
    return streams.begin();
  }
  TrackIterator get_end(){
    return streams.end();
  }

  int get_num_tracks(){
    return streams.size();
  }

  void clear();

  float get_max_x(){ return max_x; };
  float get_max_y(){ return max_y; };
  float get_max_z(){ return max_z; };

  float get_min_x(){ return min_x; };
  float get_min_y(){ return min_y; };
  float get_min_z(){ return min_z; };

  Streamline& operator [] (const int& idx) {
    return streams[idx];
  }

  //
  //  Prop/Scalar handlers
  //

  short get_num_props(){
    return num_props;
  }

  short get_num_scalars(){
    return num_scalars;
  }

  void add_prop();
  void add_scalar();
  void add_prop(const std::vector<float>& values);
  void add_scalar(const std::vector<std::vector<float> >& values);

  void wipe_header_xfm(){
    for(int i = 24; i < 36; ++i)
      header[i] = '\0';
    for(int i = 440; i < 504; ++i)
      header[i] = '\0';
  }

  //because...trackvis is dumb?
  void add_const_scalar(const std::vector<float>& values);

  void remove_prop(const int& idx);
  void remove_scalar(const int& idx);
  //void rename_prop(int& idx,std::string name);
  //void rename_scalar(int& idx,std::string name);

  void resample_trajectories(const int& num_points);


protected:
  int num_tracks;
  short num_props;
  //std::vector<std::string> prop_names;
  short num_scalars;
  //std::vector<std::string> scalar_names;

  char header[1000];

  //num_tracks x track_length (may be irregular) x 3 (xyz)
  std::vector<Streamline> streams;

  float max_x;
  float max_y;
  float max_z;

  float min_x;
  float min_y;
  float min_z;

  int verbose;

};

//class TrackIterator{
//
//};

/// I know this is just a carbon copy of std::vector now,
/// but maybe later we'll need to add our own stuff?
/// For now it'll stay simple, but it's better (IMO) than
/// baking in the std::vector representation of a streamline
///  -dmoyer 170402 1634
class Streamline{
public:
  Streamline(){;}

  void add(point pt){
    points.push_back(pt);
  }

  void delete_point(int idx){
    points.erase(points.begin() + idx);
  }

  //
  //  prop get
  //    (other prop handlers are PRIVATE, only to be called by Trackset)

  float get_prop(const int& idx){
    return props[idx];
  }

  void set_prop(const int& idx, const float& value){
    props[idx] = value;
  }

  //probably slow?
  //std::vector<float> get_scalar(
  float get_scalar(const int& idx, const int& point_idx){
    return scalars[ points.size() * idx + point_idx ];
  }

  void set_scalar(const int& idx, const int& point_idx, const float& value){
    scalars[ points.size() * idx + point_idx ] = value;
  }

  //add scalar goes here.
  //we should use flattened vector notation, for fast add/removal of scalars

  //
  //  iterators
  //

  StreamlineIterator get_iterator(){
    return points.begin();
  }

  StreamlineIterator get_end(){
    return points.end();
  }

  int size(){
    return points.size();
  }

  point& operator [] (const int& idx) {
    return points[idx];
  }

  void resample(const int& num_points);

protected:

  friend Trackset;

  //
  //    Only the Trackset class should be calling these functions.
  //    I'm sure there will be boundary cases, but if it really gets to
  //    that I'm also sure a) you know what you're doing and b) you can
  //    just modify this class.
  //
  //  prop handlers
  //

  void add_prop(const float& prop){
    props.push_back(prop);
  }

  void remove_prop(const int& idx){
    props.erase(props.begin() + idx);
  }

  void add_scalar(const std::vector<float>& values);

  void remove_scalar(const int& idx);

  void set_all_scalars(const std::vector<float>& s){
    scalars = s;
  }

  std::vector<point> points;
  std::vector<float> props;

  std::vector<float> scalars;
};


}//end of namespace trk

#endif


