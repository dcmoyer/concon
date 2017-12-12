
bin/c3_main \
  Surface_Intersect \
  --verbose 10 \
  --subj sample \
  --original_mesh_prefix "test/data/sample_subj/" \
  --original_mesh_lh_postfix "_lh_mni.m" \
  --original_mesh_rh_postfix "_rh_mni.m" \
  --trk_prefix "test/data/sample_subj/" \
  --trk_postfix "_prob.trk" \
  --output_prefix "test/data/sample_subj/outputs/" \
  --output_postfix "_xing_coords.tsv" \
  --extension_distance 4


bin/c3_main \
  bary_to_sphere \
  --verbose 10 \
  --subj sample \
  --sphere_prefix "test/data/sample_subj/" \
  --lh_sphere_postfix "_lh_mni.m" \
  --rh_sphere_postfix "_rh_mni.m" \
  --bary_coords_prefix "test/data/sample_subj/outputs/" \
  --bary_coords_postfix "_xing_coords.tsv" \
  --sphere_coords_prefix "test/data/sample_subj/outputs/" \
  --sphere_coords_postfix "_xing_sphere_coords.tsv"

bin/c3_main \
  Compute_Kernel \
  --verbose 10 \
  --subj sample \
  --sigma 0.005 \
  --epsilon 0.001 \
  --final_thold 0.000000001 \
  --OPT_VAL_exp_num_kern_samps 6 \
  --OPT_VAL_exp_num_harm_samps 5 \
  --OPT_VAL_num_harm 33 \
  --LOAD_xing_path "test/data/sample_subj/outputs/" \
  --LOAD_xing_postfix "_xing_sphere_coords.tsv" \
  --LOAD_grid_file "iso5.m" \
  --LOAD_kernel_path "test/data/sample_subj/outputs/" \
  --LOAD_kernel_postfix "_0.005.raw" \
  --LOAD_mask_file MASK \
  --SAVE_Compute_Kernel_prefix "test/data/sample_subj/outputs/"\
  --SAVE_Compute_Kernel_postfix "_0.005_full.raw" 


