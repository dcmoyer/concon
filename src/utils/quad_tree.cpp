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


#include "quad_tree.hpp"


namespace quad_tree{

//
//  Quad_Tree class
//  Front facing class to hold and deal with nodes
//

Quad_Tree::Quad_Tree(){
  max_y = max_x = min_x = min_y = -1;
  root = NULL;
}

Quad_Tree::Quad_Tree(double _max_y, double _max_x, double _min_y, double _min_x){

  max_y = _max_y; max_x = _max_x; min_y = _min_y; min_x = _min_x;

  root = new Quad_Tree_Node( NULL, _max_y, _max_x, _min_y, _min_x);

}

Quad_Tree::~Quad_Tree(){
  delete root;
}

//TODO: switch to `in_range` function
void Quad_Tree::insert(const double x, const double y, const Quad_Tree_Data data){
  assert(x >= min_x && x <= max_x && y >= min_y && y <= max_y);

  Quad_Tree_Node * current = root;
  while( ! current->is_leaf ) { //while the current node is not a leaf
    //down the rabbit hole
    if(x < current->northwest->max_x){ // if on West Side
      if(y >= current->northwest->min_y){ // North
        current = current->northwest;
      } else { //south
        current = current->southwest;
      }
    } else { // East Side
      if(y >= current->northwest->min_y){ // North
        current = current->northeast;
      } else { //south
        current = current->southeast;
      }
    }
  }

  if(current->is_empty){ //if we've found an empty leaf
    current->is_empty = false;
    current->data = data;
    return;
  }

  //otherwise it's not empty, so we must split
  Quad_Tree_Data split_data = current->data;
  do{
    //split leaf into 4
    current->is_leaf = false;
    current->northwest = new Quad_Tree_Node( current,
      current->max_y, (current->max_x + current->min_x) / 2.0,
      (current->max_y + current->min_y) / 2.0, current->min_x);
    current->southwest = new Quad_Tree_Node( current,
      (current->max_y + current->min_y)/ 2.0, (current->max_x + current->min_x) / 2.0,
      current->min_y, current->min_x);
    current->northeast = new Quad_Tree_Node( current,
      current->max_y, current->max_x,
      (current->max_y + current->min_y)/ 2.0, (current->max_x + current->min_x) / 2.0);
    current->southeast = new Quad_Tree_Node( current,
      (current->max_y + current->min_y) / 2.0, current->max_x,
      current->min_y, (current->max_x + current->min_x) / 2.0);

    //set data in the right place
    if(split_data.x < current->northwest->max_x){ // if on West Side
      if(split_data.y >= current->northwest->min_y){ // North
        current->northwest->data = split_data;
        current->northwest->is_empty = false;
      } else { //south
        current->southwest->data = split_data;
        current->southwest->is_empty = false;
      }
    } else { // East Side
      if(split_data.y >= current->northwest->min_y){ // North
        current->northeast->data = split_data;
        current->northeast->is_empty = false;
      } else { //south
        current->southeast->data = split_data;
        current->southeast->is_empty = false;
      }
    }

    //down the rabbit hole
    if(x < current->northwest->max_x){ // if on West Side
      if(y >= current->northwest->min_y){ // North
        current = current->northwest;
      } else { //south
        current = current->southwest;
      }
    } else { // East Side
      if(y >= current->northwest->min_y){ // North
        current = current->northeast;
      } else { //south
        current = current->southeast;
      }
    }
  }while(! current->is_empty); //if not empty, keep going!

  //actually add the data
  //if(current->is_empty){ //if we've found an empty leaf
  current->is_empty = false;
  current->data = data;
  //}

  return;
}

void Quad_Tree::get_range(const double _min_x, const double _min_y,
  const double _max_x, const double _max_y, std::vector<Quad_Tree_Data> &data){

  std::stack<Quad_Tree_Node *> stack;
  stack.push(root);

  Quad_Tree_Node* current;
  while(stack.size() != 0){
    current = stack.top();
    stack.pop();
    if(current->is_leaf){ //if it contains data, 
      if(!current->is_empty &&
        current->data.x < _max_x && current->data.x >= _min_x &&
        current->data.y < _max_y && current->data.y >= _min_y)
        data.push_back(current->data);
    } else { //add the children to the stack
      stack.push(current->northwest);
      stack.push(current->southwest);
      stack.push(current->northeast);
      stack.push(current->southeast);
    }

  }

  return;
}

void Quad_Tree::get_range_periodic(const double _min_x, const double _min_y,
  const double _max_x, const double _max_y, std::vector<Quad_Tree_Data> &data,
  const bool periodic_in_x, const bool periodic_in_y){

  std::vector<double> partition_x;
  std::vector<double> partition_y;

  if(periodic_in_x){
    if(_min_x < min_x){
      partition_x.push_back( max_x - (min_x - _min_x) );
      partition_x.push_back( max_x );
      partition_x.push_back( min_x );
    } else {
      partition_x.push_back( _min_x);
    }

    if(_max_x > max_x){
      partition_x.push_back( max_x );
      partition_x.push_back( min_x );
      partition_x.push_back( min_x + (_max_x - max_x));
    } else {
      partition_x.push_back( _max_x);
    }
  } else {
    partition_x = {_min_x, _max_x};
  }

  if(periodic_in_y){
    if(_min_y < min_y){
      partition_y.push_back( max_y - (min_y - _min_y) );
      partition_y.push_back( max_y );
      partition_y.push_back( min_y );
    } else {
      partition_y.push_back( _min_y);
    }

    if(_max_y > max_y){
      partition_y.push_back( max_y );
      partition_y.push_back( min_y );
      partition_y.push_back( min_y + (_max_y - max_y));
    } else {
      partition_y.push_back( _max_y);
    }
  } else {
    partition_y = {_min_y, _max_y};
  }
  int sx = (int) partition_x.size();
  int sy = (int) partition_y.size();
  for(int i = 0; i < sx; i += 2)
    for(int j = 0; j < sy; j += 2)
      get_range( partition_x[i], partition_y[j],
        partition_x[i+1], partition_y[j+1], data);

  return;
}

//
//  Quad_Tree_Node class
//  The workhorse of this thing
//

//safety first
Quad_Tree_Node::Quad_Tree_Node(){
  up = northwest = southwest = northeast = southeast = NULL;
  max_y = max_x = min_x = min_y = 0;
  is_leaf = true;
  is_empty = true;
}

Quad_Tree_Node::Quad_Tree_Node(Quad_Tree_Node * _up,
  double _max_y, double _max_x, double _min_y, double _min_x){

  up = _up;
  max_y = _max_y;
  max_x = _max_x;
  min_y = _min_y;
  min_x = _min_x;
  northwest = southwest = northeast = southeast = NULL;
  is_leaf = true;
  is_empty = true;

}

Quad_Tree_Node::~Quad_Tree_Node(){
  if(northwest != NULL) delete northwest;
  if(southwest != NULL) delete southwest;
  if(northeast != NULL) delete northeast;
  if(southeast != NULL) delete southeast;
}

void Quad_Tree_Node::set(Quad_Tree_Data _data){
  data = _data;
  is_empty = false;
}

bool Quad_Tree_Node::in_range(const double x, const double y){
  return (x >= min_x && x < max_x && y >= min_y && y < max_y);
}

bool Quad_Tree_Node::intersects_range(const double _max_y, const double _max_x,
  const double _min_y, const double _min_x){
  bool y_overlap = (min_y < _max_y && max_y >= _max_y) ||
    (max_y > _min_y && min_y <= _min_y);
  bool x_overlap = (min_x < _max_x && max_x >= _max_x) ||
    (max_x > _min_x && min_x <= _min_x);
  return y_overlap && x_overlap;
}


}//end namespace



