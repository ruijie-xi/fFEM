# fFEM

[English](README.md)

使用 Fortran 编写的二维有限元科研与教学代码，直接实现网格拓扑、有限元空间、数值积分、装配和稀疏线性代数。

## 功能

- 三角形 P0/P1/P2 元和间断线性元。
- 四边形 Q0/Q1/Q2、Nédélec 和 Raviart–Thomas 元。
- 基于内置 UMFPACK 2.0 的串行直接求解，以及 VTK 输出。
- Poisson、Stokes、多重网格和磁场示例。
- 回归测试与制造解收敛验证。

## 快速开始

需要 GNU Fortran、MPI Fortran 编译器包装器（`mpif90`）、GNU Make 和 Python 3。默认构建不需要 PETSc。

```sh
git clone https://github.com/ruijie-xi/fFEM.git
cd fFEM
make all
make check
```

`make check` 构建示例，并在临时目录运行规模受限的测试。有限元核心为串行实现；使用 `mpif90` 不代表支持分布式有限元计算。

## 目录

| 目录 | 内容 |
| --- | --- |
| `src/` | 网格、单元、装配、代数、求解器和可视化 |
| `test/` | 回归测试与收敛示例 |
| `apps/` | 多重网格课程项目和磁场应用 |
| `3rdparty/` | 内置 UMFPACK 和 Triangle 源码 |

本项目面向科研与教学，并非生产级求解器。PETSc 后端仍属实验性功能，任意网格的鲁棒性和磁场守恒性尚未充分验证。接口约定、验证范围与已知限制见[开发说明（英文）](docs/development.md)。
