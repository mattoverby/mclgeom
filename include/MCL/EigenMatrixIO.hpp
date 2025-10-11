// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_MATRIXIO_HPP
#define MCL_MATRIXIO_HPP 1

#include <Eigen/Dense>
#include <fstream>

namespace mcl {

template<typename Derived>
static inline void
write_eigen_matrix(const std::string& fn, const Eigen::MatrixBase<Derived>& A)
{
    std::ofstream file(fn);
    if (file.is_open()) {
        file << A;
    }
}

template<typename T, int dim = Eigen::Dynamic>
static inline bool
read_eigen_matrix(const std::string& fn, Eigen::Matrix<T, Eigen::Dynamic, dim>& A)
{
    // dim most often -1 (Eigen::Dynamic)

    std::vector<std::vector<T>> vals;
    std::ifstream ifs;
    ifs.open(fn.c_str());
    if (!ifs) {
        return false;
    }
    int cols = 0;

    std::string line;
    while (std::getline(ifs, line)) {
        vals.emplace_back(std::vector<T>());
        std::istringstream ss(line);
        while (ss.good()) {
            T val;
            ss >> val;
            vals.back().emplace_back(val);
        }
        cols = std::max(cols, int(vals.back().size()));
    }
    ifs.close();

    int rows = vals.size();
    if (rows == 0 || cols == 0) {
        // File exists, it's just empty
        A = Eigen::Matrix<T, Eigen::Dynamic, dim>();
        return true;
    }
    cols = dim > 0 ? std::min(cols, dim) : cols;
    A = Eigen::Matrix<T, Eigen::Dynamic, dim>::Zero(rows, cols);
    int nx = vals.size();
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < cols; ++j) {
            A(i, j) = vals[i][j];
        }
    }
    return true;

}; // end read matrix

} // end namespace mcl

#endif