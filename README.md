# Slab Reactor Toy Models

---

## Description of System

The system which I have chosen to build this toy model around is a slab reactor, which is infinite in two dimensions, and finite along the x axis, with vacuum at either end of the reactor. The reactor itself is composed of a core which is a multiplying medium, and a reflector at either side which is non-multiplying. A depiction of this system along the x axis is presented below.

![reactor](https://github.com/HunterBelanger/reactor_slab/blob/master/docs/reacteur.png)

## Installation
To make the ray-tracing and two delta-tracking programs, simply run make to build the executables in the main directory. To build the no virtual collision weighted models, you must link to the include directory, so if in this root drectory, one would run `$ clang++ -Iinclude -fopenmp nvc_delta_tracking/scatter/nvc_delta_tracking.cpp -o nvcdt`.
