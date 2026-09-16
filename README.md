# fFEM

[中文](README.zh-CN.md)

A two-dimensional finite-element research and teaching code written in Fortran.
It provides an explicit implementation of mesh topology, finite-element spaces,
quadrature, assembly, and sparse linear algebra.

## Features

- Triangular P0/P1/P2 and discontinuous linear elements.
- Quadrilateral Q0/Q1/Q2, Nédélec, and Raviart–Thomas elements.
- Serial direct solves with bundled UMFPACK 2.0 and VTK output.
- Poisson, Stokes, multigrid, and magnetic-field examples.
- Regression tests and manufactured-solution convergence checks.

## Quick start

Requires GNU Fortran, an MPI Fortran wrapper (`mpif90`), GNU Make, and Python 3.
PETSc is not required for the default build.

```sh
git clone https://github.com/ruijie-xi/fFEM.git
cd fFEM
make all
make check
```

`make check` builds the examples and runs bounded tests in temporary directories.
The core FEM implementation is serial; using `mpif90` does not imply distributed
FEM support.

## Layout

| Directory | Contents |
| --- | --- |
| `src/` | Meshes, elements, assembly, algebra, solvers, and visualization |
| `test/` | Regression tests and convergence examples |
| `apps/` | Multigrid coursework and magnetic-field applications |
| `3rdparty/` | Bundled UMFPACK and Triangle sources |

This is research software, not a production solver. The PETSc backend is
experimental; arbitrary-mesh robustness and magnetic conservation are not fully
verified. See [development notes](docs/development.md) for API contracts,
verification coverage, and limitations.
