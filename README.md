# GeoPH: Delaunay–Rips persistence diagrams

This repository contains the code for the manuscript [*Delaunay–Rips filtration: a study and an algorithm*](https://arxiv.org/abs/2512.17382).

It also contains the code for the planar Rips and Delaunay--Rips persistence algorithm described in [*Topological Autoencoders++: Fast and Accurate Cycle-Aware Dimensionality Reduction*](https://arxiv.org/abs/2502.20215).

## Dependencies

The following dependencies are required:
* Boost >= 1.72
* CGAL
* Eigen3

For multi-threading, the following dependencies are additionally required:
* Boost >= 1.83
* OpenMP
* TBB

This implementation also uses and contains the following codes:
* [ankerl::unordered_dense](https://github.com/martinus/unordered_dense), a densely stored hashmap and hashset by Martin Leitner-Ankerl
* [dset](https://github.com/wjakob/dset), a lock-free parallel disjoint set data structure by Wenzel Jakob

## C++

This code can be used as a header-only C++ library. Copy and include the `geoPHUtils.h` and `geoPH{2,3,d}.h` files (in the `src` directory) in your project. Then use the classes `DRPersistence{2,3,D}`.

## Python

Execute
`pip wheel .` in the project root
and
`pip install` the produced `.whl` wheel file.
