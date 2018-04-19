# concon

This is an implementation of *Con*tinuous *Con*nectivity from
[1](#references)
and
[2](#references).
It is meant to be a parcellation free representation of
brain connectivity, i.e. a form of connecitivty outside of the discrete
network frame. If none of that made any sense to you, we strongly recommend
you, we strongly recommend you read the papers and email the first author
(Daniel Moyer, moyerd [at] usc.edu) with any questions.  

This project will be expanding as time goes on, but for now is clearly in
pre-alpha! Sorry! DM will add tutorials and hopefully real documentation
soon.

# Requirements

-A gcc compiler that supports c++11   
-cmake   
-OpenMP  
-The Eigen library from [http://eigen.tuxfamily.org/](http://eigen.tuxfamily.org/). This is only used by third party software.

# Installation

Step 0: Install required software from [Requirements](#requirements). Don't
forget to add Eigen's library to your `CPLUS_INCLUDE_PATH`, or equivalent.
This is usually done on Unix systems by adding the following to your .bashrc:

    export CPLUS_INCLUDE_PATH="/usr/include/eigen3/:${CPLUS_INCLUDE_PATH}"

For non-Ubuntu systems this path may vary.  

Step 1: Create a build directory wherever you would like to build concon.  

Step 2: Navigate to the build directory and run
`cmake /path/to/concon/repo/`, with the appropriate path.  

Step 3: Run `make`.

The binaries and static libraries for concon will be built in `bin` and `lib`.

# Third Party Software

We include the spherical harmonics package from
[here](https://github.com/google/spherical-harmonics). Their 
original build instructions used Bazel; we have created a short cmake recipe.
It is otherwise unchanged. Please see their License file for the corresponding
copyright information.

# License

Please see our License file (MIT).

# References

[1] D Moyer et al. *Continuous Representations of Brain Connectivity using Spatial Point Processes*, Medical Image Analysis (MedIA) 2017.  

[2] D Moyer et al. *A Continuous Model of Cortical Connectivity*, International Conference on Medical Image Computing and Computer-Assisted Intervention (MICCAI) pp. 157-165. Springer International Publishing, 2016.  

