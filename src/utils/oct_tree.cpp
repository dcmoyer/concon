
#include "oct_tree.hpp"


namespace oct_tree {


void Oct_Tree::init(const std::pair<double,double>& x_min_max,
    const std::pair<double,double>& y_min_max,
    const std::pair<double,double>& z_min_max,
    const double _min_bin_size){

  init_flag = true;
  x_min = x_min_max.first;
  x_max = x_min_max.second;
  y_min = y_min_max.first;
  y_max = y_min_max.second;
  z_min = z_min_max.first;
  z_max = z_min_max.second;

  min_bin_size = _min_bin_size;

  root = new OT_Node();
  root->set(NULL,
    x_min, x_max,
    y_min, y_max,
    z_min, z_max);
}

void Oct_Tree::init(const double _x_min, const double _x_max,
    const double _y_min, const double _y_max,
    const double _z_min, const double _z_max,
    const double _min_bin_size){

  init_flag = true;
  x_min = _x_min;
  x_max = _x_max;
  y_min = _y_min;
  y_max = _y_max;
  z_min = _z_min;
  z_max = _z_max;

  min_bin_size = _min_bin_size;

  root = new OT_Node();
  root->set(NULL,
    x_min, x_max,
    y_min, y_max,
    z_min, z_max);

}

Oct_Tree::~Oct_Tree(){

  if(root != NULL)
    delete root;

}

int helper_q(bool x_flag, bool y_flag, bool z_flag){

  if(z_flag){
    //1,2,3,4
    if(y_flag){
      //1,2
      if(x_flag){
        return 1;
      } else {
        return 2;
      }
    } else {
      //3,4
      if(x_flag){
        return 4;
      } else {
        return 3;
  }}} else {
    //5,6,7,8
    if(y_flag){
      //5,6
      if(x_flag){
        return 5;
      } else {
        return 6;
      }
    } else {
      //7,8
      if(x_flag){
        return 8;
      } else {
        return 7;
  }}}
}


void Oct_Tree::insert(const double x, const double y, const double z,
    int idx){

  //TODO range check

  if(!init_flag){
    std::cerr << "[oct_tree] "
      << "insert called without init, aborting." << std::endl;
    exit(1);
  }

  OT_Node * target = root;
  bool x_flag, y_flag, z_flag;
  while(! (target->is_leaf) ){
    x_flag = x >= (target->x_max + target->x_min)/2.0;
    y_flag = y >= (target->y_max + target->y_min)/2.0;
    z_flag = z >= (target->z_max + target->z_min)/2.0;

#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
    target = &(target->q[helper_q(x_flag, y_flag, z_flag) - 1]);
    //target = target->q[helper_q(x_flag, y_flag, z_flag) - 1];
  }

  bool add_flag = ((target->x_max - target->x_min)/2.0 < min_bin_size ||
      (target->y_max - target->y_min)/2.0 < min_bin_size ||
      (target->z_max - target->z_min)/2.0 < min_bin_size) &&
        !(target->is_empty);
  //push along the current items
  while(! (target->is_empty)){
    //TODO: check for exact matches...
    if(
      (target->x_max - target->x_min)/2.0 < min_bin_size ||
      (target->y_max - target->y_min)/2.0 < min_bin_size ||
      (target->z_max - target->z_min)/2.0 < min_bin_size ){
      add_flag = true;
      break;
    }
    x_flag = x >= (target->x_max + target->x_min)/2.0;
    y_flag = y >= (target->y_max + target->y_min)/2.0;
    z_flag = z >= (target->z_max + target->z_min)/2.0;
    //target = target->q + (helper_q(x_flag, y_flag, z_flag) - 1);
    target->split();
    target = &(target->q[helper_q(x_flag, y_flag, z_flag) - 1]);
  }

  OT_Data d;
  d.x = x;
  d.y = y;
  d.z = z;
  d.idx = idx;
  d.next = NULL;

  if(add_flag){
    target->add_to_list(d);
  } else {
    target->set(d);
  }
}

//get
void Oct_Tree::get_range(const double _x_min, const double _x_max,
  const double _y_min, const double _y_max,
  const double _z_min, const double _z_max,
  std::vector<Oct_Tree::OT_Data>& results){

  if(!init_flag){
    std::cerr << "[oct_tree] "
      << "get_range called without init, aborting." << std::endl;
    exit(1);
  }

  get_range_helper( root,
    _x_min, _x_max,
    _y_min, _y_max,
    _z_min, _z_max,
    results);

  return;
}

void Oct_Tree::get_range_helper(OT_Node * target,
  const double _x_min, const double _x_max,
  const double _y_min, const double _y_max,
  const double _z_min, const double _z_max,
  std::vector<Oct_Tree::OT_Data>& results){

  if(target == NULL){
    return;
  }

  if(!target->intersects_range(
    _x_min, _x_max,
    _y_min, _y_max,
    _z_min, _z_max)){
    return;
  }

  if(!target->is_leaf && !target->is_empty){
    std::cout << "WAT A NON_LEAF IS FILLED" << std::endl;
    exit(1);
  }

  if(! target->is_empty){
    if(target->data_in_range(
      _x_min, _x_max,
      _y_min, _y_max,
      _z_min, _z_max)){

      //std::cerr << "Line 203, something is wrong here." << std::endl;
      //exit(1);

      if(target->is_list){
        //TODO: SOMETHING IS WRONG HERE
        //modified for lists
        Oct_Tree::OT_Data * next_list_obj = target->data;
        while(next_list_obj != NULL){
          results.push_back(* next_list_obj );
          next_list_obj = next_list_obj->next;
        }
      } else {
        results.push_back(* (target->data) );
      }
    }
  }
  if(target->is_leaf){
    return;
  }

  for(int i = 0; i < 8; ++i){
    get_range_helper(target->q+i,
      _x_min, _x_max,
      _y_min, _y_max,
      _z_min, _z_max,
      results);
  }

}

void Oct_Tree::print(){
  std::list<std::pair<OT_Node*,int> > the_queue;
  the_queue.push_back(std::pair<OT_Node*,int>(root,0));

  if(!init_flag){
    std::cout << "Oct Tree Is Empty" << std::endl;
    return;
  }

  std::pair<OT_Node*,int> p;
  OT_Node* current;
  while(!the_queue.empty()){

    p = the_queue.front();
    current = p.first;
    the_queue.pop_front();

    std::cout << "Level " << p.second <<
      " x:[" << current->x_min << "," << current->x_max << "] " <<
      " y:[" << current->y_min << "," << current->y_max << "] " <<
      " z:[" << current->z_min << "," << current->z_max << "] " <<
      std::endl;

    if(!current->is_empty){
      std::cout << "\tData: idx:" << current->data->idx 
        << " x:" << current->data->x
        << " y:" << current->data->y
        << " z:" << current->data->z
        << std::endl;
    }

    if(!current->is_leaf){
      for(int i = 0; i < 8; ++i){
        the_queue.push_back(
          std::pair<OT_Node*,int>(current->q + i,p.second + 1)
        );
      }
    }
  }
}



//
// OT_Node class
//

OT_Node::OT_Node(){
  is_empty = true;
  is_leaf = true;
  is_list = false;
  parent = NULL;
  q = NULL;
}

OT_Node::OT_Node(OT_Node * up, const double _x_min, const double _x_max,
    const double _y_min, const double _y_max,
    const double _z_min, const double _z_max){

  set(up,
    _x_min,_x_max,
    _y_min,_y_max,
    _z_min,_z_max
  );

}

void OT_Node::set(OT_Node * up, const double _x_min, const double _x_max,
    const double _y_min, const double _y_max,
    const double _z_min, const double _z_max){

  x_min = _x_min;
  x_max = _x_max;
  y_min = _y_min;
  y_max = _y_max;
  z_min = _z_min;
  z_max = _z_max;

  is_empty = true;
  is_leaf = true;
  parent = up;
  q = NULL;
}


OT_Node::~OT_Node(){

  if(!is_leaf){
    delete[] q;
  }

  if(!is_empty){
    clear_data();
  }
}

void OT_Node::set(Oct_Tree::OT_Data _data){
  if(!is_empty && !is_list){
    std::cout << "Trying to fill already full OT Node. Aborting." <<
      std::endl;  
    exit(1);
  }

  if(is_list){
    std::cout << "Cannot set listed node, aborting." << std::endl;
    exit(1);
  }

  data = new Oct_Tree::OT_Data;
  (*data) = _data;
  data->next = NULL;
  is_empty = false;
}

void OT_Node::add_to_list(Oct_Tree::OT_Data _data){
  if(is_empty){
    std::cout << "Trying to listify an empty node, aborting."
      << std::endl;
    exit(1);
  }

  is_list = true;

  Oct_Tree::OT_Data* target = data;
  while(target->next != NULL){
    target = target->next;
  }
  target->next = new Oct_Tree::OT_Data;
  (*(target->next)) = _data;
  target->next->next = NULL;
}

void OT_Node::clear_data(){
  if(!is_empty){
    Oct_Tree::OT_Data* target = data;
    while(target->next != NULL){
      data = target;
      target = target->next;
      delete data;
    }
  }
  data = NULL;
  is_empty = true;

}

bool helper_in_range(
  const double& t, const double& t_min, const double& t_max){
  return (t >= t_min && t < t_max);
}

bool helper_in_interval(
  const double& lower, const double& upper,
  const double& t_min, const double& t_max){
  return helper_in_range(lower, t_min, t_max) ||
    helper_in_range(upper, t_min, t_max) ||
    (lower <= t_min && t_max <= upper);
}

bool OT_Node::in_range(
  const double& x, const double& y, const double& z){
  return (
    helper_in_range(x,x_min,x_max) && 
    helper_in_range(y,y_min,y_max) && 
    helper_in_range(z,z_min,z_max)
  );
}

bool OT_Node::data_in_range(
  const double& _x_min, const double& _x_max,
  const double& _y_min, const double& _y_max,
  const double& _z_min, const double& _z_max){
  return (
    helper_in_range(data->x,_x_min,_x_max) && 
    helper_in_range(data->y,_y_min,_y_max) && 
    helper_in_range(data->z,_z_min,_z_max)
  );
}

bool OT_Node::intersects_range(
  const double& _x_min, const double& _x_max,
  const double& _y_min, const double& _y_max,
  const double& _z_min, const double& _z_max
  ){
  return (
    helper_in_interval(_x_min,_x_max,x_min,x_max) && 
    helper_in_interval(_y_min,_y_max,y_min,y_max) && 
    helper_in_interval(_z_min,_z_max,z_min,z_max)
  );

}

void OT_Node::split(){
  if( q != NULL ){
    std::cout << "attempt to split OT_Node while already split, aborting."
      << std::endl;
    exit(1);
  }

  if( is_empty ){
    std::cout << "attempt to split non-empty OT_Node, aborting" << std::endl;
    exit(1);
  }

  if( is_list ){
    std::cout << "attempt to split listed OT_Node, aborting" << std::endl;
    exit(1);
  }

  double half_x = (x_max + x_min)/2.0;
  double half_y = (y_max + y_min)/2.0;
  double half_z = (z_max + z_min)/2.0;

  q = new OT_Node[8];
  q[0] = OT_Node(this, half_x, x_max, half_y, y_max, half_z, z_max);
  q[1] = OT_Node(this, x_min, half_x, half_y, y_max, half_z, z_max);
  q[2] = OT_Node(this, x_min, half_x, y_min, half_y, half_z, z_max);
  q[3] = OT_Node(this, half_x, x_max, y_min, half_y, half_z, z_max);

  q[4] = OT_Node(this, half_x, x_max, half_y, y_max, z_min, half_z);
  q[5] = OT_Node(this, x_min, half_x, half_y, y_max, z_min, half_z);
  q[6] = OT_Node(this, x_min, half_x, y_min, half_y, z_min, half_z);
  q[7] = OT_Node(this, half_x, x_max, y_min, half_y, z_min, half_z);


  bool x_flag = data->x >= half_x;
  bool y_flag = data->y >= half_y;
  bool z_flag = data->z >= half_z;

  q[helper_q(x_flag, y_flag, z_flag) - 1].set(*data);

  clear_data();
  is_leaf = false;

}

}//end of namespace oct_tree



