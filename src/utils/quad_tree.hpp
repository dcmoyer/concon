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

