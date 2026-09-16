# Development and numerical contracts

[Repository overview](../README.md)

The default, verified backend is serial assembly with the bundled UMFPACK 2.0
direct solver.

## Build and verify

Requirements: GNU Fortran, an MPI Fortran compiler wrapper (`mpif90`, used by the
build and MPI smoke example), GNU Make, and Python 3 for the test runner. PETSc is
not required by the default build. No Python third-party packages are required.

```sh
make all
make check
```

`make check` builds the core, test programs, homework applications, and magnetic
examples. It runs bounded algebra/finite-element regressions and manufactured
solution convergence checks. Numerical output is written into a temporary
directory, not into tracked experimental data. Example executables can also be
built using `make -C apps/magnetic all` or `make -C apps/homework/final all`.
When running examples directly, create their requested output directories first.

Makefile builds are deliberately serialized. Core module dependencies are
conservative: changing a core source rebuilds the core objects. Run `make
clean-all` before changing compiler/optimization flags or switching PETSc on/off.

## API and numerical contracts

- `VECTOR` stores `size` and allocatable `data`; a fresh vector has size zero.
  `VectorInit`, `VectorReset`, and `VectorCopy` accept fresh destinations.
  Legacy `v%Init(n)`, `v%Reset(n)`, `v%Norm()`, and `v%AddVector(w[,a])` remain
  available. `Norm()` and `VectorNormLinf` are infinity norms; the historical
  `VectorNormL2` is an RMS coefficient norm, not an FE integral norm. Use
  `ComputeNorm`/`ComputeError` for finite-element integral quantities.
- `Interpolate` and FE assembly/evaluation accept `type(VECTOR)`, not raw arrays.
  `ComputeDof` is a function. Scalar derivatives of vector fields, such as RT
  divergence, return one component in all three FE evaluation entry points.
- Triplet `N_nz` is capacity; `actual_nnz` is the number of assembled entries.
  Duplicate entries are summed. Addition uses actual entries only. Conversion
  to CSC removes exact zeros, **not** small nonzero coefficients. Numerical
  dropping is opt-in: `call MatrixColumnTrim(A, tolerance)`.
- Reference triangle integration weights sum to `1/2`; quadrilateral weights
  sum to `1`. Vertex-based rules are integration rules, not normalized averages.
- `call SolverSolveUMFPACK2(A,b,x[,status])` accepts triplet or CSC storage and
  initializes/resizes `x`. Inputs must be a square nonempty matrix and a matching
  RHS. Zero status means success, including valid duplicate-triplet assembly.
  Singular matrices and other factorization/solve failures produce nonzero
  status and invalidate `x` with NaNs. Without `status`, failure is an `error stop`.
  Input and output vectors must be distinct actual arguments.
- `AddBdryMarker` accepts an optional tolerance (default `1d-8`).

## Verification coverage

The regression suite checks small-entry preservation, rectangular transpose,
triplet spare capacity, vector API compatibility, solver residuals and failure
paths, diagonal/block-triangular solves, triangle quadrature, default/explicit
boundary tolerances, constant-field reproduction for ten spaces (including
nonaffine quadrilaterals), and RT point/edge divergence.

Convergence checks use three mesh levels:

| Problem | Required behavior |
| --- | --- |
| Q2 Poisson | approximately order 3 in L2, order 2 in H1 |
| P2/P1 Stokes | approximately order 3/2 for velocity L2/H1; order 2 for pressure L2 |
| Quadrilateral RT1 interpolation | order 1 in L2 and H(div) |
| Quadrilateral RT2 interpolation | order 2 in L2 and H(div) |

RT1/RT2 are this repository's names for the four-/twelve-local-DOF spaces.
They are not a claim about a universal polynomial-degree naming convention.
Stokes uses a pressure gauge rather than an artificial pressure mass term.
Multigrid Jacobi/Gauss-Seidel examples are checked for residual reduction.
Magnetic moving/static examples receive bounded smoke tests, **not** a complete
conservation or temporal-convergence certification. `test_mag` is compile-checked;
its external mesh/input-dependent scenario is not run by the default suite.

## Scope and remaining limitations

These tests establish a regression baseline, not mathematical correctness for
all meshes, coefficients, boundary conditions, or applications. Supply valid,
nondegenerate, consistently oriented meshes. Arbitrary/inverted mesh validation,
full element polynomial/commuting-diagram certification, general boundary
constraint handling, and magnetic conservation tests remain future work.

The optional PETSc implementation is experimental and not covered by `make check`.
In particular, its distributed ownership and lifecycle handling have not been
repaired or validated here; do not treat it as a supported MPI FEM solver.
The default Stokes example now uses the verified UMFPACK path. UMFPACK workspace
allocation remains conservative and factorizations are not cached; this patch
does not address large-scale performance or add 3D/adaptive capabilities.

## Bundled solver correction

The active split UMFPACK source and its monolithic source both initialize
`iout`/`xout` in `ums2f0`. Singleton-only block-triangular factorizations can skip
the routines that otherwise assign these flags; reading them uninitialized
caused false workspace failures. This is a localized initialization correction,
not a change to the factorization algorithm. The bundled archive now uses
replacement (`ar rcs`) instead of accumulating duplicate object members.
