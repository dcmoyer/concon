# concon

This is an implementation of continuous connectivity from [1] and [2].

[1] D Moyer et al. *Continuous Representations of Brain Connectivity using Spatial Point Processes*, Medical Image Analysis (MedIA) 2017.  

[2] D Moyer et al. *A Continuous Model of Cortical Connectivity*, International Conference on Medical Image Computing and Computer-Assisted Intervention (MICCAI) pp. 157-165. Springer International Publishing, 2016.  


Requirements
----
-A gcc compiler that supports c++11   
-OpenMP  
-The Eigen library from [http://eigen.tuxfamily.org/](http://eigen.tuxfamily.org/). This is only used by third party software.

Third Party Software
----

We include the spherical harmonics package from
[here](https://github.com/google/spherical-harmonics). Their 
original build instructions used Bazel; we have created a short cmake recipe.
It is otherwise unchanged. Please see their License file for the corresponding
copyright information (Apache 2).

License
----

Please see our License file (MIT).


