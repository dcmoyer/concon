
#ifndef TRK_HPP
#define TRK_HPP

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <functional>

//TODO:remove this in favor of homebrew'd one
#include <string.h>

namespace trk {

const int TRK_LOCAL_ARRAY_SIZE = 1500000;

class Tractset {
public:
  Tractset();
  Tractset(int _verbose);
  ~Tractset();

  void set_verbose(int _verbose){
    verbose = _verbose;
  }

  //int load_ucf(std::string filename);
  int load_trk(std::string filename);

  int save_trk(std::string filename);

  //access
  int get_num_tracts(){
    return num_tracts;
  }

  int get_num_points(const int& idx){
    if(idx > num_tracts){
      std::cerr << "Bad idx lookup in trk. Aborting." << std::endl;
      exit(1);
    }
    return num_points[idx];
  }

  void clear();

  //dangerous right here
  //bad practice, but maybe we need to go FAST. idk, takes out an extra
  //function call every time we need to access a tract.
  void get_references(float*** & _tracts_x_coords, int* & _num_points ){
    _tracts_x_coords = tracts_x_coords;
    _num_points = num_points;
  }

  void set_scalar(const int& trk_idx, const int& pt_idx, const int& s_idx,
    const float& value);
  void set_prop(const int& trk_idx, const int& p_idx, const float& value);

  //this one's a dangerous one, initialized to zero though.
  void add_scalar();
  void add_prop();

  short get_num_scalars(){
    return num_scalars;
  }

  short get_num_props(){
    return num_props;
  }

  float get_scalar(const int& trk_idx, const int& pt_idx, const int& s_idx){
    return tracts_x_scalars[trk_idx][pt_idx][s_idx];
  }
  float get_prop(const int& trk_idx, const int& p_idx){
    return tracts_x_props[trk_idx][p_idx];
  }

  //maybe this should be a C style array? To force people to be good people?
  //Morality is hard tho.
  //void add_scalar(std::vector<std::vector<float>> scalar);

  void add_prop(const std::vector<float>& prop);

  //This function is because TrackVis can't get it's format to work with its
  //viewer.
  void add_trk_const_scalar(std::vector<float> scalar);

  void remove_scalar(const int& idx);
  void remove_prop(const int& idx);

  //why is this here? Because it needs direct access to the tracts, and
  //because I thought it would be better than writing it in python.
  //Perhaps I should have written it in Python.
  //But oh well.
  //
  //This function resamples trajectories to N points, even if N is
  //larger than the number of points in the original tract. This is for use
  //with MDP distance.
  //
  //TODO?:Split off into their own trk child classes?
  void resample_trajectories(int N);

  void resample_tracts(
    int N,
    bool replacement=false,
    bool in_order=true,
    int seed=1919
  );

  //really really inefficient?
  //this will leave a trailing list of pointers to NOTHING
  //at the end of it, which will only dealloc on clear() or exit.
  void delete_tract(int idx);

  //TRUE if KEPT, KEPT if TRUE
  //THIS DELETES THINGS!
  void filter_by_const_scalar(int prop_idx,
    std::function<bool (float)> condition);
  void filter_by_prop(int prop_idx,
    std::function<bool (float)> condition);

  //TRUE if KEPT, KEPT if TRUE
  //THIS DELETES THINGS!
  //void filter_by_scalar(bool (*condition)(float* values, int len));

  //since I don't have code to do streamline tracing, this will be left null
  // for a while
  //void add_tract();

//if this is too slow, make public and use direct access
protected:

  //trk only stores single precision
  float*** tracts_x_coords;

  //scalars held along the tracts
  short num_scalars;
  float*** tracts_x_scalars;

  //scalars help per tract
  short num_props;
  float** tracts_x_props;

  //vector of vectors might be too big?

  char header[1000];

  int num_tracts;
  int* num_points;
  int verbose;

};


}//end of namespace trk

#endif

