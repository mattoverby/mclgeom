# mclgeom

Header-only functions and tools for nonlinear optimization, physics-based animation, and general mesh processing.
Most of this code was implemented throughout my PhD research, so please use the citations within the file if used. Documentation is a WIP.

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

Quasi-Newton optimizer for general nonlinear minimization problems.

```cpp
mcl::LBFGS<MatrixXd> solver;
// One function must be defined for both the objective and/or gradient.
// The gradient computation can be skipped if x.size() != g.size().
solver.gradient = [&](const MatrixXd& x, MatrixXd& g) -> double {
    double obj = compute_objective(x);
    if (g.rows() == x.rows()) g = compute_gradient(x);
    return obj;
};
double result = solver.minimize(x);
```

### Karush-Kuhn-Tucker (KKT) Solver

Conjugate gradient-based solver for saddle-point systems (linear constrained optimization).

```cpp
mcl::KKTSolver<VectorXd, SparseMatrixXd> kkt;
int iters = kkt.solve(A, b, C, d, x, y);
```

### Multi-Color Gauss-Seidel

Parallel iterative solver for `Ax = b` with optional projection operators. Uses graph coloring for thread-safe updates.

```cpp
mcl::MultiColorGaussSeidel<MatrixXd> mcgs;
mcgs.project = [](int i, MatrixXd& X) { /* optionally enforce constraints */ };
mcl::graph_color(A, colors);
int iters = mcgs.solve(A, B, X, colors);
```

### Levenberg-Marquardt

Damped least-squares solver for underdetermined systems `min ||f(x)||^2`. Suitable for over-determined residual problems.

```cpp
mcl::LevenbergMarquardt<VectorXd, SparseMatrixXd> lm;
lm.objective = [&](const VectorXd& x, VectorXd& residual, SparseMatrixXd& J, bool need_J) {
    /* compute residual and optionally Jacobian J */
};
double result = lm.iterate(x);
```

## Physics-Based Animation

### Energy Models

Nonlinear material energy densities for deformable bodies. Includes neo-Hookean, St. Venant-Kirchhoff, ARAP, and symmetric Dirichlet.

```cpp
using Model = mcl::StableNeoHookean<3, double>;
double energy = Model::energy_density(lame_params, F);
VectorXd gradient; // call Model::gradient(lame_params, F, gradient)
MatrixXd hessian; // call Model::hessian(lame_params, F, hessian)
```

### Bending Models

Cloth and surface bending energies. Includes quadratic Bergou model (Q) and linear Volino-Magnenat-Thalmann model (K0).

```cpp
mcl::make_hinges(F, H); // Extract hinge edges from triangle mesh
auto Q = mcl::quadratic_bend_Q(x0, x1, x2, x3);
auto K0 = mcl::quadratic_bend_K0(x0, x1, x2, x3);
```

## License

**mclgeom** MIT license.

