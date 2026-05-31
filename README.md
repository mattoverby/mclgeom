# mclgeom

Header-only functions and tools for nonlinear optimization, physics-based animation, and general mesh processing.
Most of this code was implemented throughout my PhD research, so please use the citations within if used.

By Matt Overby ([https://mattoverby.net](https://mattoverby.net))

## Build

![Tests](https://github.com/mattoverby/mclgeom/actions/workflows/build_and_test.yml/badge.svg)

Nearly every file is self-contained and can simply be dropped into your own programs. Many of the components use the following dependencies:
- Eigen ([https://gitlab.com/libeigen/eigen](https://gitlab.com/libeigen/eigen))
- Thread Building Blocks ([https://github.com/uxlfoundation/oneTBB](https://github.com/uxlfoundation/oneTBB))

To build all tests and GUI examples:
```sh
mkdir build
cd build
cmake .. -DMCLGEOM_EXAMPLES=ON -DMCLGEOM_UNIT_TESTS=ON -DCMAKE_BUILD_TYPE=Release
make -j
./test/test_lbfgs
```

## Numerical Optimization

There are several general purpose solvers for common problems in computer graphics.

### L-BFGS

### Karush-Kuhn-Tucker (KKT) Matrix Solver

### Multi-Color Gauss-Seidel

### Levenberg-Marquardt

## Physics-Based Animation

### Energy Density Functions

### Cloth Bending Models

### Impact Zones

## License

**mclgeom** MIT license.

