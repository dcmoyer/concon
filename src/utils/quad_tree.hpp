//
//
//  quad_tree.hpp
//
//  dmoyer 160215 1525
//
//
//  conventions:
//    North/South == +/- Y-Dim
//    East/West == +/- X-Dim
//    Data \in [min_*, max_*) HALF OPEN interval
//

#ifndef QUAD_TREE_HPP
#define QUAD_TREE_HPP

#include <iostream>
#include <vector>
#include <cassert>
#include <stack>

namespace quad_tree{

class Quad_Tree_Node;

struct Quad_Tree_Data {
  int idx;
  double x;
  double y;
};

class Quad_Tree {
public:
  Quad_Tree();
  Quad_Tree(double _max_y, double _max_x, double _min_y, double _min_x);
  ~Quad_Tree();

  void insert(const double x, const double y, const Quad_Tree_Data data);
  
  //int get(const double x, const double y);
  void get_range(const double _min_x, const double _min_y,
    const double _max_x, const double _max_y, std::vector<Quad_Tree_Data> &data);

  //this is a special function for periodic domains
  void get_range_periodic(const double _min_x, const double _min_y,
    const double _max_x, const double _max_y, std::vector<Quad_Tree_Data> &data,
    const bool periodic_in_x=true, const bool periodic_in_y=true);

private:
  double max_y;
  double max_x;
  double min_y;
  double min_x;

  Quad_Tree_Node * root;

  friend Quad_Tree_Node;

};


class Quad_Tree_Node{
public:
  Quad_Tree_Node();
  Quad_Tree_Node(Quad_Tree_Node * _up,
    double _max_y, double _max_x, double _min_y, double _min_x);
  ~Quad_Tree_Node();

  void set(Quad_Tree_Data _data);

  bool in_range(const double x, const double y);

  bool intersects_range(const double _max_y, const double _max_x,
    const double _min_y, const double _min_x);


private:
  Quad_Tree_Node * up;
  Quad_Tree_Node * northwest;
  Quad_Tree_Node * southwest;
  Quad_Tree_Node * northeast;
  Quad_Tree_Node * southeast;

  double max_y;
  double max_x;
  double min_y;
  double min_x;

  //int data;
  //TODO: redo as pointer to conserve space
  Quad_Tree_Data data;
  bool is_leaf;
  bool is_empty;

  friend Quad_Tree;

};




}//end namespace
#endif

