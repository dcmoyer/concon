//
//
// oct_tree.hpp
//
// dmoyer 170327 1040
//
// Why do I reimplement things? Because I'm a sad grad student, that's why.
//
//conventions:
//  q1: x > 0, y > 0, z > 0
//  q2: x < 0, y > 0, z > 0
//  q3: x < 0, y < 0, z > 0
//  q4: x > 0, y < 0, z > 0
//
//     y+
//  q2 | q1
//  ---z--- x+
//  q3 | q4
//
//  q5: x > 0, y > 0, z > 0
//  q6: x < 0, y > 0, z > 0
//  q7: x < 0, y < 0, z > 0
//  q8: x > 0, y < 0, z > 0
//
//     y+
//  q6 | q5
//  ---z--- x+
//  q7 | q8
//

#ifndef OCT_TREE_HPP
#define OCT_TREE_HPP

#include <iostream>
#include <vector>
#include <cassert>
#include <stack>
#include <list>
#include <utility>

namespace oct_tree {

class OT_Node;

class Oct_Tree {
public:
  struct OT_Data {
    int idx;
    double x;
    double y;
    double z;
    OT_Data* next;
  };

  Oct_Tree(){
    init_flag = false;
    min_bin_size = -1;
  }

  //STL why not
  void init(const std::pair<double,double>& x_min_max,
    const std::pair<double,double>& y_min_max,
    const std::pair<double,double>& z_min_max,
    const double _min_bin_size=-1);

  void init(const double _x_min, const double _x_max,
    const double _y_min, const double _y_max,
    const double _z_min, const double _z_max,
    const double _min_bin_size=-1);

  ~Oct_Tree();

  void insert(const double x, const double y, const double z,
    int idx);

  void get_range(const double _x_min, const double _x_max,
    const double _y_min, const double _y_max,
    const double _z_min, const double _z_max,
    std::vector<Oct_Tree::OT_Data>& results);

  void print();

private:

  void get_range_helper(OT_Node * target,
    const double _x_min, const double _x_max,
    const double _y_min, const double _y_max,
    const double _z_min, const double _z_max,
    std::vector<Oct_Tree::OT_Data>& results);

  bool init_flag;
  double min_bin_size;

  double x_min;
  double x_max;
  double y_min;
  double y_max;
  double z_min;
  double z_max;

  OT_Node * root;

  friend OT_Node;

};

class OT_Node {
public:
  OT_Node();
  OT_Node(OT_Node * up, const double _x_min, const double _x_max,
    const double _y_min, const double _y_max,
    const double _z_min, const double _z_max);
  void set(OT_Node * up, const double _x_min, const double _x_max,
    const double _y_min, const double _y_max,
    const double _z_min, const double _z_max);
  ~OT_Node();

  void set(Oct_Tree::OT_Data _data);
  void add_to_list(Oct_Tree::OT_Data _data);

  void clear_data();

  bool in_range(const double& x, const double& y, const double& z);

  bool intersects_range(
    const double& _x_min, const double& _x_max,
    const double& _y_min, const double& _y_max,
    const double& _z_min, const double& _z_max
  );

  bool data_in_range(
    const double& _x_min, const double& _x_max,
    const double& _y_min, const double& _y_max,
    const double& _z_min, const double& _z_max
  );

  void split();

private:
  double x_min;
  double x_max;
  double y_min;
  double y_max;
  double z_min;
  double z_max;
  OT_Node * q;

  Oct_Tree::OT_Data * data;
  OT_Node * parent;
  bool is_empty;
  bool is_leaf;
  bool is_list;

  friend Oct_Tree;

};

}//end of namespace oct_tree

#endif

