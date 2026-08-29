// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_GEOM_CONSTRAINT_POLISH_HPP
#define MCL_GEOM_CONSTRAINT_POLISH_HPP 1

#include <MCL/ConstraintZone.hpp>
#include <MCL/LevenbergMarquardt.hpp>

#include <Eigen/Core>

#include <iostream>
#include <set>
#include <unordered_set>

namespace mcl {

/// @brief For solving sets of inequality constraints
/// @tparam T scalar type
/// @tparam DIM dimension of vertices, with primitive dimension DIM + 1
template<typename T, int DIM>
class ConstraintSolver
{
  public:
    std::unordered_set<int> primitives_in_set; ///< primitive indices of volume constraints
    std::vector<VolumeConstraint<T, DIM>> volume_constraints;
    std::vector<ConstraintZone> zones; ///< constraints that share a vertex
    bool needs_merge = true;           ///< if constraint set has changed

    /// @brief Clears existing data
    void clear()
    {
        primitives_in_set.clear();
        volume_constraints.clear();
        zones.clear();
        needs_merge = true;
    }

    /// @brief Returns number of constraints
    size_t num_constraints() const { return volume_constraints.size(); }

    /// @brief Loops over all primitives and adds any that are inverted to the set
    void add_inversions(const T* x, const T* x_rest, const int* primitives, int num_primitives)
    {
        T threshold = T(1e-4);
        volume_constraints.reserve(num_primitives / 4);
        primitives_in_set.reserve(num_primitives / 4);
        zones.reserve(num_primitives / 4);
        for (int i = 0; i < num_primitives; ++i) {

            // If the primitive is already in the constraint set, no need to re-check
            if (primitives_in_set.count(i) > 0) {
                continue;
            }

            // Any new constraint is initially added to its own zone
            auto stencil = get_primitive<DIM + 1>(i, primitives);
            auto verts = get_verts<T, DIM, DIM + 1>(x, stencil.data());
            int constraint_index = volume_constraints.size();
            if constexpr (DIM == 2) {
                if (signed_triangle_area(verts[0], verts[1], verts[2]) <= threshold) {
                    needs_merge = true;
                    volume_constraints.emplace_back(x_rest, stencil, -1);
                    primitives_in_set.emplace(i);
                    zones.emplace_back(constraint_index, stencil.data(), stencil.size());
                }
            } else if constexpr (DIM == 3) {
                if (signed_tet_volume(verts[0], verts[1], verts[2], verts[3]) <= threshold) {
                    needs_merge = true;
                    volume_constraints.emplace_back(x_rest, stencil, -1);
                    primitives_in_set.emplace(i);
                    zones.emplace_back(constraint_index, stencil.data(), stencil.size());
                }
            }
        }
    }

    /// @brief Moves vertices to best satisfy all constraints.
    /// @return Number of iterations and updates vertices (x)
    int solve(T* x, const T* x_rest, int num_vertices, const int* primitives, int num_primitives)
    {
        // Initial gather of constraints
        add_inversions(x, x_rest, primitives, num_primitives);
        if (volume_constraints.empty()) {
            return 0;
        }

        // Merge constraints into zones
        if (needs_merge) {
            needs_merge = false;
            ConstraintZone::merge_zones(num_vertices, zones);
        }

        // Solver loop
        int iter = 0;
        int max_iterations = 20;
        while (iter < max_iterations) {
            ++iter;

            // Solve constraints
            for (auto& zone : zones) {
                iterate_zone(zone, x);
            }

            // Check for new constraints
            add_inversions(x, x_rest, primitives, num_primitives);

            // Check for termination
            if (check_termination(x)) {
                break;
            }

            // Combine new constraints into new zones
            if (needs_merge) {
                needs_merge = false;
                ConstraintZone::merge_zones(num_vertices, zones);
            }
        }

        return iter;
    }

    /// @brief Computes a (local) delta x to minimize constraint residuals
    void iterate_zone(ConstraintZone& zone, T* global_x)
    {
        using VectorType = Eigen::VectorX<T>;
        using MatrixType = Eigen::SparseMatrix<T>;
        LevenbergMarquardt<VectorType, MatrixType> LM;

        LM.objective = [&](const VectorType& local_x, VectorType& r, MatrixType& J, bool needJ) -> void {
            std::vector<Eigen::Triplet<T>> J_triplets;
            if (needJ) {
                J_triplets.reserve(zone.constraints.size() * (DIM + 1) * DIM);
                J.resize(zone.constraints.size(), local_x.rows());
                J.setZero();
            }

            int r_index = -1;
            r = VectorType::Zero(zone.constraints.size());
            for (auto constraint_index : zone.constraints) {
                r_index++;

                // Remap global stencil to local stencil
                const auto& constraint = volume_constraints[constraint_index];
                std::array<int, DIM + 1> local_stencil;
                for (int i = 0; i < DIM + 1; ++i) {
                    local_stencil[i] = zone.global_to_local[constraint.stencil[i]];
                }

                const auto verts = get_verts<T, DIM, DIM + 1>(local_x.data(), local_stencil.data());
                T eval = constraint.eval(verts);
                if (eval >= constraint.target_eval()) {
                    continue;
                }

                if (needJ) {
                    const auto gradients = constraint.gradients(verts);
                    for (int i = 0; i < DIM + 1; ++i) {
                        for (int j = 0; j < DIM; ++j) {
                            assert(local_stencil[i] * DIM + j < local_x.rows());
                            J_triplets.emplace_back(r_index, local_stencil[i] * DIM + j, gradients[i][j]);
                        }
                    }
                }

                r[r_index] = eval - constraint.target_eval();
            }

            if (!J_triplets.empty()) {
                J.setFromTriplets(J_triplets.begin(), J_triplets.end());
            }
        };

        // Get local vertices
        VectorType local_x = VectorType::Zero(zone.stencil.size() * DIM);
        for (size_t i = 0; i < zone.stencil.size(); ++i) {
            int global_index = zone.stencil[i];
            for (int j = 0; j < DIM; ++j) {
                local_x[i * DIM + j] = global_x[global_index * DIM + j];
            }
        }

        // Solve
        LM.iterate(local_x);
        for (size_t i = 0; i < zone.stencil.size(); ++i) {
            int global_index = zone.stencil[i];
            for (int j = 0; j < DIM; ++j) {
                global_x[global_index * DIM + j] = local_x[i * DIM + j];
            }
        }
    }

    bool check_termination(const T* global_x)
    {
        for (auto& c : volume_constraints) {
            const auto verts = get_verts<T, DIM, DIM + 1>(global_x, c.stencil.data());
            if (c.eval(verts) < std::min(T(0), c.target_eval())) {
                return false;
            }
        }
        return true;
    }
};

} // end ns mcl

#endif // MCL_GEOM_CONSTRAINT_POLISH_HPP