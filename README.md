# mclgeom

Header-only functions and tools for nonlinear optimization, physics-based animation, and general mesh processing.
Most of this code was implemented throughout my PhD research.

By Matt Overby ([https://mattoverby.net](https://mattoverby.net))

## Build

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

## Nonlinear Optimization



## Physics-Based Animation

### Energy Density Functions

## Mesh Processing


