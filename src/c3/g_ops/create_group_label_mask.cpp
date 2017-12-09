
#include "c3/g_ops/create_group_label_mask.hpp"


namespace c3 {

int create_group_label_mask( std::vector<std::string> subj_list,
  std::unordered_map<std::string,std::string> params){

  std::string label_prefix, label_lh_postfix, label_rh_postfix;
  std::string grid_prefix, grid_lh_postfix, grid_rh_postfix;
  std::string out_prefix, out_postfix;
  bool verbose = false;
  bool and_flag = false;
  bool face_mask_flag = false;
  std::vector<int> target_regions;
  int thold = -1;

  for( std::pair<std::string,std::string> p : params ){
    if(p.first == "verbose"){
      verbose = true;
    }

    if(p.first == "group_label_regions"){
      std::string temp = p.second;
      std::stringstream ss;
      ss << temp;
      while(!ss.eof()){
        int temp_target;
        ss >> temp_target;
        target_regions.push_back(temp_target);
        ss.get();

        if(verbose)
          std::cout << temp_target << " ";
      }

      if(verbose)
        std::cout << std::endl;
    }

    if(p.first == "label_prefix"){
      label_prefix = p.second;
    }

    if(p.first == "label_lh_postfix"){
      label_lh_postfix = p.second;
    }

    if(p.first == "label_rh_postfix"){
      label_rh_postfix = p.second;
    }

    if(p.first == "grid_prefix"){
      grid_prefix = p.second;
    }

    if(p.first == "grid_lh_postfix"){
      grid_lh_postfix = p.second;
    }

    if(p.first == "grid_rh_postfix"){
      grid_rh_postfix = p.second;
    }

    if(p.first == "group_label_out_prefix"){
      out_prefix = p.second;
    }

    if(p.first == "group_label_out_postfix"){
      out_postfix = p.second;
    }

    if(p.first == "group_label_AND"){
      and_flag = true;
    }

    if(p.first == "face_mask_flag"){
      face_mask_flag = true;
    }

    if(p.first == "label_mask_thold"){
      thold = std::stoi(p.second);
    }
  }

  std::string region_string;
  std::string and_or_string;
  for( int& region : target_regions ){

    if(region == -1){
      region_string = "medial_wall";
    } else {
      region_string = std::to_string(region);
    }

    if(and_flag){
      and_or_string = "AND";
    } else {
      and_or_string = "OR";
    }
    std::vector<bool> labels;
    std::vector<int> label_counts;

    bool first = true;
    for( std::string &subj : subj_list ){

      if(labels.size() == 0){
        first = true;
      } else {
        first = false;
      }

      std::string lh_file = label_prefix + subj + label_lh_postfix;
      std::ifstream input_lh_file(lh_file.c_str());
      std::string rh_file = label_prefix + subj + label_rh_postfix;
      std::ifstream input_rh_file(rh_file.c_str());
      int temp;
      int max = -1;
      int idx = 0;
      bool check = false;

      if(!input_rh_file.good() || !input_lh_file.good()){
        continue;
      }

      while(input_rh_file >> temp){
        if(temp > max)
          max = temp;
        if(first){
          labels.push_back(temp == region);
          label_counts.push_back(temp == region);
        } else if(and_flag){
          labels[idx] = labels[idx] && (temp == region);
          label_counts[idx] += (int) (temp == region);
        } else {
          labels[idx] = labels[idx] || (temp == region);
          label_counts[idx] += (int) (temp == region);
        }
        idx += 1;
        check = check || (temp == region);
      }
      input_rh_file.close();

      while(input_lh_file >> temp){
        if(temp >= 0){
          temp = temp + max;
        }

        if(first){
          labels.push_back(temp == region);
          label_counts.push_back(temp == region);
        } else if(and_flag){
          labels[idx] = labels[idx] && (temp == region);
          label_counts[idx] += (int) (temp == region);
        } else {
          labels[idx] = labels[idx] || (temp == region);
          label_counts[idx] += (int) (temp == region);
        }
        idx += 1;
        check = check || (temp == region);
      }
      input_lh_file.close();

      if(!check){
        std::cout << subj << std::endl;
      }
    }

    //
    //  Face Mask
    //

    std::string face_string = "";
    std::vector<int> face_mask;
    if(face_mask_flag){

      face_string = "FACE" + std::to_string(thold);

      MeshLib::Solid rh_grid, lh_grid;
      std::string grid_file = grid_prefix + grid_rh_postfix;
      rh_grid.read( grid_file.c_str() );
      grid_file = grid_prefix + grid_lh_postfix;
      lh_grid.read( grid_file.c_str() );
      //masking for the medial wall
      for(int grid_idx = 0; grid_idx < 2; ++grid_idx){
        //...I guess there's a more elegant way
        //but today is not that day
        int diff = 0;
        int face_diff = 0;

        MeshLib::Solid * target_grid;
        if(grid_idx == 0){
          target_grid = &rh_grid;
        } else {
          target_grid = &lh_grid;
          diff = rh_grid.numVertices();
          face_diff = rh_grid.numFaces();
        }

        for(MeshLib::SolidFaceIterator face_iter(target_grid); !face_iter.end();
          ++face_iter){

          //check surrounding faces
          bool surrounded = true;
          for(MeshLib::FaceVertexIterator vertex_iter(face_iter.value());
            !vertex_iter.end(); ++vertex_iter){
            surrounded = surrounded &&
              (label_counts[vertex_iter.value()->id() + diff - 1] >= thold);
          }

          //remember the 1 indexing ffs
          if(surrounded){
            face_mask.push_back(thold);
          } else {
            face_mask.push_back(0);
          }
        }
      }

      //face mask to vertex mask
      label_counts = face_mask;
    }

    std::ofstream mask_output( out_prefix + face_string +
      "group_" + region_string + "_" +
      and_or_string + out_postfix);

    //for( bool label_value : labels ){ 
    for( int label_value : label_counts ){ 
      mask_output << label_value << std::endl;
    }
    mask_output.close();
  }

  return 0;
}





}// end of namespace c3


