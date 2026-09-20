// Copyright Matt Overby 2026.
// Distributed under the MIT License.

#ifndef MCL_GEOM_MESH_INJECTIVITY_SOLVER_HPP
#define MCL_GEOM_MESH_INJECTIVITY_SOLVER_HPP

#include "ConstraintZone.hpp"
#include "LevenbergMarquardt.hpp"

#include <Eigen/Core>
#include <tbb/parallel_for.h>

#include <iostream>
#include <unordered_set>

namespace mcl {

/// @brief A robust solver for global injectivity constraints.
/// Used to compute foldover-free maps for mesh parameterization and deformation.
/// See Overby et al. 2021 (https://doi.org/10.1111/cgf.14361) for details.
/// This is a reimplementation from the original code release and I haven't fully vetted it yet.
/// TODO: Collision constraints.
/// @tparam T scalar type
/// @tparam DIM dimension of vertices, with primitive dimension DIM + 1
template<typename T, int DIM>
class MeshInjectivitySolver
{
  protected:
    std::unordered_set<int> primitives_in_set; ///< primitive indices of volume constraints
    std::unordered_set<int> pinned_vertices;   ///< Dirichlet boundary condition
    std::vector<VolumeConstraint<T, DIM>> volume_constraints;
    std::vector<ConstraintZone> zones; ///< constraints that share a vertex
    bool needs_merge = true;           ///< if constraint set has changed
    int max_iterations = 1000;

  public:
    /// @brief Clears existing data
    void clear()
    {
        primitives_in_set.clear();
        pinned_vertices.clear();
        volume_constraints.clear();
        zones.clear();
        needs_merge = true;
    }

    /// @brief Sets Dirichlet boundary conditions on the vertices
    void add_pins(const int* pinned_vertex_inds, int num_pinned_vertices)
    {
        for (int i = 0; i < num_pinned_vertices; ++i) {
            pinned_vertices.emplace(pinned_vertex_inds[i]);
        }
    }

    /// @brief Moves vertices to best satisfy all constraints. Note the solver attempts to enforce the target
    /// volume for all elements, but will exit once all elements have a positive volume.
    /// @return Number of iterations and updates vertices (x)
    int solve(T* x,
              const T* x_rest, // may be null if rest unavailable
              int num_vertices,
              const int* primitives,
              int num_primitives)
    {
        // Initial gather of constraints, one zone per constraint
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
        while (iter < max_iterations) {
            ++iter;

            // Solve constraints in parallel.
            tbb::parallel_for(tbb::blocked_range<int>(0, int(zones.size())), [&](const tbb::blocked_range<int>& range) {
                for (int i = range.begin(); i != range.end(); ++i) {
                    iterate_zone(zones[i], x);
                }
            });

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

  private:
    /// @brief Returns true if the vertex is fixed/pinned.
    bool is_fixed_vertex(int index) { return pinned_vertices.count(index) > 0; }

    /// @brief Returns true if a stencil has at least one non-pinned vertex
    bool has_free_vertex(const Eigen::Vector<int, DIM + 1>& stencil)
    {
        for (int i = 0; i < DIM + 1; ++i) {
            if (!is_fixed_vertex(stencil[i])) {
                return true;
            }
        }
        return false;
    }

    /// @brief Loops over all primitives and adds any that are inverted to the set. x_rest may be null.
    void add_inversions(const T* x, const T* x_rest, const int* primitives, int num_primitives)
    {
        T threshold = T(1e-8); // target_eval?
        volume_constraints.reserve(num_primitives / 4);
        primitives_in_set.reserve(num_primitives / 4);
        zones.reserve(num_primitives / 4);
        for (int i = 0; i < num_primitives; ++i) {

            // If the primitive is already in the constraint set, no need to re-check
            if (primitives_in_set.count(i) > 0) {
                continue;
            }

            auto stencil = get_primitive<DIM + 1>(i, primitives);
            if (!has_free_vertex(stencil)) {
                continue;
            }

            // Any new constraint is initially added to its own zone
            auto verts = get_verts<T, DIM, DIM + 1>(x, stencil.data());
            int constraint_index = volume_constraints.size();
            if constexpr (DIM == 2) {
                if (signed_triangle_area(verts[0], verts[1], verts[2]) <= threshold) {
                    needs_merge = true;
                    volume_constraints.emplace_back(x_rest == nullptr ? x : x_rest, stencil, -1);
                    primitives_in_set.emplace(i);
                    zones.emplace_back(constraint_index, stencil.data(), stencil.size());
                }
            } else if constexpr (DIM == 3) {
                if (signed_tet_volume(verts[0], verts[1], verts[2], verts[3]) <= threshold) {
                    needs_merge = true;
                    volume_constraints.emplace_back(x_rest == nullptr ? x : x_rest, stencil, -1);
                    primitives_in_set.emplace(i);
                    zones.emplace_back(constraint_index, stencil.data(), stencil.size());
                }
            }
        }
    }

    /// @brief Computes a (local) delta x to minimize constraint residuals
    void iterate_zone(ConstraintZone& zone, T* global_x)
    {
        using VectorType = Eigen::VectorX<T>;
        using MatrixType = Eigen::SparseMatrix<T>;

        // TODO: Keep around LM parameter for each zone
        LevenbergMarquardt<VectorType, MatrixType> LM;

        // Reuse J_triplets/r_values buffer to avoid repeated allocation
        std::vector<Eigen::Triplet<T>> J_triplets;
        J_triplets.reserve(zone.constraints.size() * (DIM + 1) * DIM);
        std::vector<T> r_values;
        r_values.reserve(zone.constraints.size());

        // Get local vertices, removing fixed vertices. Otherwise, we end up with
        // zeros on the matrix diagonal. This is a little convoluted because now we have two
        // different local-to-global mappings: one for the zone (which contains fixed vertices)
        // and one for the LM solver (no fixed vertices).
        // Might be more efficient to cache these maps/reduced spaces and only update
        // when the zone has changed.
        std::vector<int> LM_local_to_global; // no fixed vertices in this map
        std::unordered_map<int, int> LM_global_to_local;
        std::vector<T> local_x_data;
        LM_local_to_global.reserve(zone.stencil.size());
        local_x_data.reserve(zone.stencil.size() * DIM);
        for (size_t i = 0; i < zone.stencil.size(); ++i) {
            int global_index = zone.stencil[i];
            if (!is_fixed_vertex(global_index)) {
                LM_global_to_local.emplace(global_index, LM_local_to_global.size());
                LM_local_to_global.emplace_back(global_index);
                for (int j = 0; j < DIM; ++j) {
                    local_x_data.emplace_back(global_x[global_index * DIM + j]);
                }
            }
        }

        // Function to compute the objective and/or Jacobian.
        LM.objective = [&](const VectorType& local_x, VectorType& r, MatrixType& J, bool needJ) -> void {
            J_triplets.clear();
            r_values.clear();

            for (auto constraint_index : zone.constraints) {

                const auto& constraint = volume_constraints[constraint_index];

                // The LM solver's variable, local_x, does NOT contain fixed variables
                // and uses a different local-global mapping than the constraint zone.
                std::array<Eigen::Vector<T, DIM>, DIM + 1> verts;
                for (int i = 0; i < DIM + 1; ++i) {
                    int global_index = constraint.stencil[i];
                    if (is_fixed_vertex(global_index)) {
                        verts[i] = Eigen::Map<const Eigen::Vector<T, DIM>>(global_x + global_index * DIM).eval();
                    } else {
                        int local_index = LM_global_to_local[global_index];
                        verts[i] = local_x.template segment<DIM>(local_index * DIM).eval();
                    }
                }

                // Only include the active set in the residual and Jacobian
                T eval = constraint.eval(verts);
                if (eval >= constraint.target_eval()) {
                    continue; // r = 0
                }

                if (needJ) {
                    int J_row = r_values.size();
                    auto gradients = constraint.gradients(verts);
                    for (int i = 0; i < DIM + 1; ++i) {
                        int global_index = constraint.stencil[i];
                        if (!is_fixed_vertex(global_index)) {
                            int local_index = LM_global_to_local[global_index];
                            for (int j = 0; j < DIM; ++j) {
                                J_triplets.emplace_back(J_row, local_index * DIM + j, gradients[i][j]);
                            }
                        }
                    }
                }

                r_values.emplace_back(eval - constraint.target_eval());
            }

            if (!r_values.empty()) {
                r = Eigen::Map<VectorType>(r_values.data(), r_values.size()).eval();
                if (!J_triplets.empty()) {
                    J.resize(r_values.size(), local_x.size());
                    J.setFromTriplets(J_triplets.begin(), J_triplets.end());
                }
            }
        };

        VectorType local_x = Eigen::Map<VectorType>(local_x_data.data(), local_x_data.size()).eval();

        // Solve
        T objective = LM.iterate(local_x);

        // Map back to global buffer if there wasn't an error.
        if (objective >= 0) {
            for (size_t i = 0; i < LM_local_to_global.size(); ++i) {
                int global_index = LM_local_to_global[i];
                for (int j = 0; j < DIM; ++j) {
                    global_x[global_index * DIM + j] = local_x[i * DIM + j];
                }
            }
        }
    }

    /// @brief Checks if all constraints are sufficiently solved.
    bool check_termination(const T* global_x, bool target_positive_volume = true)
    {
        for (auto& c : volume_constraints) {
            T target_volume = target_positive_volume ? std::numeric_limits<T>::epsilon() : c.target_eval();
            const auto verts = get_verts<T, DIM, DIM + 1>(global_x, c.stencil.data());
            if (c.eval(verts) < target_volume) {
                return false;
            }
        }
        return true;
    }
};

} // end ns mcl

#endif // MCL_GEOM_MESH_INJECTIVITY_SOLVER_HPP