// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_GEOM_CONSTRAINT_POLISH_HPP
#define MCL_GEOM_CONSTRAINT_POLISH_HPP 1

#include "SignedMeasure.hpp"
#include "DisjointSets.hpp"
#include "LevenbergMarquardt.hpp"

#include <Eigen/Core>

#include <unordered_set>
#include <iostream>

namespace mcl {

/// @brief A volume constraint, volume(x) >= target
/// @tparam T scalar type
/// @tparam DIM dimension of vertices, with primitive dimension DIM + 1
template<typename T, int DIM>
class VolumeConstraint
{
public:
    Eigen::Vector<int, DIM+1> stencil = Eigen::Vector<int, DIM+1>::Zero();
    int index = -1; ///< unique index set by solver
    T scaling = 1;

    /// @brief Constructor
    VolumeConstraint() = default;

    /// @brief Constructor, scaling is computed as the surface area if not positive
    VolumeConstraint(const T* x, const Eigen::Vector<int, DIM+1> stencil_, int index_, T scaling_ = 1);

    /// @brief Destructor
    ~VolumeConstraint() = default;

    /// @brief Returns (minimum) target volume
    T target_eval() const { return 1e-3; }

    /// @brief Returns nonlinear eval at x
    T eval(const std::array<Eigen::Vector<T,DIM>, DIM+1> &verts) const;

    /// @brief Returns gradients at x
    std::array<Eigen::Vector<T,DIM>, DIM+1> gradients(const std::array<Eigen::Vector<T,DIM>, DIM+1> &verts) const;
};

/// @brief For containing and organizing inequality constraints
/// @tparam T scalar type
/// @tparam DIM dimension of vertices, with primitive dimension DIM + 1
template<typename T, int DIM>
class InequalityConstraintSet
{
public:
    struct ConstraintZone
    {
        int index = -1; ///< unique index of this zone
        std::vector<int> constraints; ///< constraint index
        std::vector<int> vertex_mapping; ///< local -> global vertex index
    };

    std::unordered_set<int> primitives_in_set; ///< primitive indices of volume constraints
    std::vector<VolumeConstraint<T,DIM>> volume_constraints;
    std::vector<ConstraintZone> zones; ///< constraints that share a vertex
    std::vector<std::pair<int,int>> global_to_local; ///< vertex -> [zone, local vertex]
    bool needs_sort = true; ///< if constraint set has changed

    /// @brief Clears existing data
    void clear()
    {
        primitives_in_set.clear();
        volume_constraints.clear();
        zones.clear();
        global_to_local.clear();
        needs_sort = true;
    }

    /// @brief Returns number of constraints
    size_t num_constraints() const
    {
        return volume_constraints.size();
    }

    /// @brief Loops over all primitives and adds any that are inverted to the set
    void add_inversions(const T* x, const T* x_rest, const int* primitives, int num_primitives)
    {
        T threshold = T(1e-4);
        needs_sort = true;
        int constraint_index = volume_constraints.size();
        volume_constraints.reserve(num_primitives/4);
        primitives_in_set.reserve(num_primitives/4);
        for (int i=0; i<num_primitives; ++i) {

            // If the primitive is already in the constraint set, no need to re-check
            if (primitives_in_set.count(i) > 0) {
                continue;
            }

            auto stencil = get_primitive<DIM+1>(i, primitives);
            auto verts = get_verts<T,DIM,DIM+1>(x, stencil.data());
            if constexpr (DIM == 2) {
                if (signed_triangle_area(verts[0], verts[1], verts[2]) <= threshold) {
                    volume_constraints.emplace_back(x_rest, stencil, constraint_index++, -1);
                    primitives_in_set.emplace(i);
                }
            }
            else if constexpr (DIM == 3) {
                if (signed_tet_volume(verts[0], verts[1], verts[2], verts[3]) <= threshold) {
                    volume_constraints.emplace_back(x_rest, stencil, constraint_index++, -1);
                    primitives_in_set.emplace(i);
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
        sort(num_vertices);

        // Solver loop
        int iter = 0;
        int max_iterations = 20;
        while(iter < max_iterations) {
            ++iter;

            // Solve constraints
            for (auto &zone : zones) {
                iterate_zone(zone, x);
            }

            // Check for new constraints
            add_inversions(x, x_rest, primitives, num_primitives);

            // Check for termination
            if (check_termination(x)) {
                break;
            }

            // Combine new constraints into new zones
            sort(num_vertices);
        }

        return iter;
    }

    /// @brief Group all constraints that share vertices, i.e., constraint (impact) zone.
    void sort(int num_vertices)
    {
        // TODO: Retain zones and merge new ones into existing zones
        // to retain LM parameter?

        needs_sort = false;
        zones.clear();
        global_to_local.clear();
        if (num_constraints() == 0) {
            return;
        }

        // Disjoint set to find connected stencils
        int max_index = -1;
        DisjointSets dj(num_vertices);
        for (auto &c : volume_constraints) {
            max_index = std::max(max_index, c.stencil[0]);
            for (int i = 1; i<DIM+1; ++i) {
                dj.make_union(c.stencil[i], c.stencil[i-1]);
                max_index = std::max(max_index, c.stencil[i]);
            }
        }

        // Sort into zones
        std::vector<std::unordered_set<int>> zone_stencils;
        std::vector<int> parent_to_zone_index;
        parent_to_zone_index.reserve(num_constraints()/4);
        for (auto &c : volume_constraints) {
            int parent = dj.find(c.stencil[0]);
            while(parent >= parent_to_zone_index.size()) {
                parent_to_zone_index.emplace_back(-1);
            }

            // new zone
            if (parent_to_zone_index[parent] < 0) {
                int new_zone_index = zones.size();
                parent_to_zone_index[parent] = new_zone_index;
                zones.emplace_back();
                zone_stencils.emplace_back();
                zones.back().index = new_zone_index;
            }

            int zone_index = parent_to_zone_index[parent];
            zones[zone_index].constraints.emplace_back(c.index);
            for (int i = 0; i<DIM+1; ++i) {
                zone_stencils[zone_index].emplace(c.stencil[i]);
            }
        }

        // Map local/global vertex indices for zones
        global_to_local.resize(num_vertices, {-1, -1});
        for (auto &zone : zones) {
            zone.vertex_mapping.resize(zone_stencils[zone.index].size());
            int local_index = 0;
            for (int global_index : zone_stencils[zone.index]) {
                zone.vertex_mapping[local_index] = global_index;
                assert(global_to_local[global_index].first == -1);
                global_to_local[global_index] = {zone.index, local_index};
                ++local_index;
            }
        }
    }

    /// @brief Computes a (local) delta x to minimize constraint residuals
    void iterate_zone(ConstraintZone &zone, T *global_x)
    {
        using VectorType = Eigen::VectorX<T>;
        using MatrixType = Eigen::SparseMatrix<T>;
        LevenbergMarquardt<VectorType, MatrixType> LM;

        LM.objective = [&](const VectorType &local_x, VectorType &r, MatrixType &J, bool needJ) -> void
        {
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
                const auto &constraint = volume_constraints[constraint_index];
                std::array<int, DIM + 1> local_stencil;
                for (int i=0; i<DIM + 1; ++i) {
                    local_stencil[i] = global_to_local[constraint.stencil[i]].second;
                    assert(global_to_local[constraint.stencil[i]].first == zone.index);
                }
                 
                const auto verts = get_verts<T,DIM,DIM+1>(local_x.data(), local_stencil.data());
                T eval = constraint.eval(verts);
                if (eval >= constraint.target_eval()) {
                    continue;
                }

                if (needJ) {
                    const auto gradients = constraint.gradients(verts);
                    for (int i=0; i<DIM + 1; ++i) {
                        for (int j=0; j<DIM; ++j) {
                            assert(local_stencil[i]*DIM+j < local_x.rows());
                            J_triplets.emplace_back(r_index, local_stencil[i]*DIM+j, gradients[i][j]);
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
        VectorType local_x = VectorType::Zero(zone.vertex_mapping.size() * DIM);
        for (size_t i =0; i < zone.vertex_mapping.size(); ++i) {
            int global_index = zone.vertex_mapping[i];
            for (int j=0; j<DIM; ++j) {
                local_x[i*DIM+j] = global_x[global_index*DIM+j];
            }
        }

        // Solve
        LM.iterate(local_x);
        for (size_t i = 0; i < zone.vertex_mapping.size(); ++i) {
            int global_index = zone.vertex_mapping[i];
            for (int j=0; j<DIM; ++j) {
                global_x[global_index*DIM+j] = local_x[i*DIM+j];
            }
        }
    }

    bool check_termination(const T *global_x) {
        for (auto &c : volume_constraints) {
            const auto verts = get_verts<T,DIM,DIM+1>(global_x, c.stencil.data());
            if(c.eval(verts) < std::min(T(0), c.target_eval())) {
                return false;
            }
        }
        return true;
    }
};

//
// Implementation
//

template<typename T, int DIM>
VolumeConstraint<T,DIM>::VolumeConstraint(const T* x, const Eigen::Vector<int, DIM+1> stencil_, int index_, T scaling_) : stencil(stencil_), index(index_), scaling(scaling_)
{
    if (scaling <= T(0)) {
        auto verts = get_verts<T,DIM,DIM+1>(x, stencil.data());
        if constexpr (DIM == 2) {
            scaling = T(1) / (T(0.5) * triangle_perimeter(verts[0], verts[1], verts[2]));
        }
        else if constexpr (DIM == 3) {
            scaling =  T(1) / (T(0.5) * tet_surface_area(verts[0], verts[1], verts[2], verts[3]));
        }
    }
}

template<typename T, int DIM>
T VolumeConstraint<T,DIM>::eval(const std::array<Eigen::Vector<T,DIM>, DIM+1> &verts) const
{
    if constexpr (DIM == 2)
    {
        return scaling * signed_triangle_area(verts[0], verts[1], verts[2]);
    }
    else if constexpr (DIM == 3)
    {
        return scaling * signed_tet_volume(verts[0], verts[1], verts[2], verts[3]);
    }
    return 0;  
}

template<typename T, int DIM>
std::array<Eigen::Vector<T,DIM>, DIM+1> VolumeConstraint<T,DIM>::gradients(const std::array<Eigen::Vector<T,DIM>, DIM+1> &verts) const
{
    if constexpr (DIM == 2)
    {
        auto grads = signed_triangle_area_gradients(verts[0], verts[1], verts[2]);
        grads[0] *= scaling;
        grads[1] *= scaling;
        grads[2] *= scaling;
        return grads;
    }
    else if constexpr (DIM == 3)
    {
        auto grads = signed_tet_volume_gradients(verts[0], verts[1], verts[2], verts[3]);
        grads[0] *= scaling;
        grads[1] *= scaling;
        grads[2] *= scaling;
        grads[3] *= scaling;
        return grads;
    }
    return {};
}


} // end ns mcl

#endif // MCL_GEOM_CONSTRAINT_POLISH_HPP