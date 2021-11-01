# Stochastic-model-of-T-cell-repositioning

# Introduction

This is the software implementing the models used by the Group of Prof. Heiko Rieger at Saarland University to study the T Cell repolarization during target elimination.
The T Cell establishes contact with one target cell and creates only one IS.

We publish it here to foster good scientific practices, i.e. to facilitate reproduction of our work by interested parties.

For more information can be found [here](https://www.rieger.uni-saarland.de/homepage/research/biological_physics/research_publications/T_Cell_modelling.html
).


# Building

* Eigen 3
* Open MPI
* g++ > 4

## Successfully built on

* 4.1.7-hardened-r1-ARCH
* 5.8.14-arch1-1

## Build process

CMake is our build system. You most likely need to configure the build according to your system specifics. Basic buils can be done by placing both folders in the build directory.

\$ cd build_directory <br />
\$ cmake one_IS_source_code/ -DCMAKE_BUILD_TYPE="" . <br />
\$ make <br />


## Details

### MPI

The MPI is NOT required and the numerical results can be collected(more slowly) without it. When simulating only one process to visualize the repositioning, MPI is unnecessary.

## Visualization

Simulation produces output files into the folder picturesVideos/textFiles. To make a movie 


\$ cd picturesVideos <br />
\$ sh Make_Movie.sh <br />




