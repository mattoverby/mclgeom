// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_GEOM_CONSTRAINT_ZONE_HPP
#define MCL_GEOM_CONSTRAINT_ZONE_HPP 1

#include "DisjointSets.hpp"
#include "SignedMeasure.hpp"

#include <Eigen/Core>

#include <set>
#include <unordered_map>

namespace mcl {

/// @brief A volume constraint, volume(x) >= target
/// @tparam T scalar type
/// @tparam DIM dimension of vertices, with primitive dimension DIM + 1
template<typename T, int DIM>
class VolumeConstraint
{
  public:
    Eigen::Vector<int, DIM + 1> stencil = Eigen::Vector<int, DIM + 1>::Zero();
    T scaling = 1;

    /// @brief Constructor
    VolumeConstraint() = default;

    /// @brief Constructor, scaling is computed as the surface area if not positive
    VolumeConstraint(const T* x, const Eigen::Vector<int, DIM + 1> stencil_, T scaling_ = 1);

    /// @brief Destructor
    ~VolumeConstraint() = default;

    /// @brief Returns (minimum) target volume
    T target_eval() const { return 1e-3; }

    /// @brief Returns nonlinear eval at x
    T eval(const std::array<Eigen::Vector<T, DIM>, DIM + 1>& verts) const;

    /// @brief Returns gradients at x
    std::array<Eigen::Vector<T, DIM>, DIM + 1> gradients(const std::array<Eigen::Vector<T, DIM>, DIM + 1>& verts) const;
};

/// @brief A collection of constraints that share a vertex
class ConstraintZone
{
  public:
    int index = -1;                               ///< unique index of this zone
    std::vector<int> constraints;                 ///< global constraint index
    std::vector<int> stencil;                     ///< local -> global vertex indices
    std::unordered_map<int, int> global_to_local; ///< global -> local vertex indices

    /// @brief Destructor
    ~ConstraintZone() = default;

    /// @brief Constructor
    ConstraintZone() = default;

    /// @brief Constructor
    ConstraintZone(int constraint_index, const int* sten, int stencil_size);

    /// @brief Merges another zone into this one
    void merge(const ConstraintZone& zone);

    /// @brief Merges all zones that share indices
    static void merge_zones(int num_vertices, std::vector<ConstraintZone>& zones);
};

//
// Implementation
//

template<typename T, int DIM>
VolumeConstraint<T, DIM>::VolumeConstraint(const T* x, const Eigen::Vector<int, DIM + 1> stencil_, T scaling_)
    : stencil(stencil_)
    , scaling(scaling_)
{
    if (scaling <= T(0)) {
        auto verts = get_verts<T, DIM, DIM + 1>(x, stencil.data());
        if constexpr (DIM == 2) {
            scaling = T(1) / (T(0.5) * triangle_perimeter(verts[0], verts[1], verts[2]));
        } else if constexpr (DIM == 3) {
            scaling = T(1) / (T(0.5) * tet_surface_area(verts[0], verts[1], verts[2], verts[3]));
        }
    }
}

template<typename T, int DIM>
T
VolumeConstraint<T, DIM>::eval(const std::array<Eigen::Vector<T, DIM>, DIM + 1>& verts) const
{
    if constexpr (DIM == 2) {
        return scaling * signed_triangle_area(verts[0], verts[1], verts[2]);
    } else if constexpr (DIM == 3) {
        return scaling * signed_tet_volume(verts[0], verts[1], verts[2], verts[3]);
    }
    return 0;
}

template<typename T, int DIM>
std::array<Eigen::Vector<T, DIM>, DIM + 1>
VolumeConstraint<T, DIM>::gradients(const std::array<Eigen::Vector<T, DIM>, DIM + 1>& verts) const
{
    if constexpr (DIM == 2) {
        auto grads = signed_triangle_area_gradients(verts[0], verts[1], verts[2]);
        grads[0] *= scaling;
        grads[1] *= scaling;
        grads[2] *= scaling;
        return grads;
    } else if constexpr (DIM == 3) {
        auto grads = signed_tet_volume_gradients(verts[0], verts[1], verts[2], verts[3]);
        grads[0] *= scaling;
        grads[1] *= scaling;
        grads[2] *= scaling;
        grads[3] *= scaling;
        return grads;
    }
    return {};
}

ConstraintZone::ConstraintZone(int constraint_index, const int* sten, int stencil_size)
{
    stencil.insert(stencil.end(), sten, sten + stencil_size);
    constraints.emplace_back(constraint_index);
    for (size_t i = 0; i < stencil.size(); ++i) {
        global_to_local[stencil[i]] = i;
    }
}

void
ConstraintZone::merge(const ConstraintZone& zone)
{
    std::set<int> combined_constraints(constraints.begin(), constraints.end());
    std::set<int> combined_stencil(stencil.begin(), stencil.end());
    combined_constraints.insert(zone.constraints.begin(), zone.constraints.end());
    combined_stencil.insert(zone.stencil.begin(), zone.stencil.end());
    constraints.assign(combined_constraints.begin(), combined_constraints.end());
    stencil.assign(combined_stencil.begin(), combined_stencil.end());
    global_to_local.clear();
    for (size_t i = 0; i < stencil.size(); ++i) {
        global_to_local[stencil[i]] = i;
    }
}

void
ConstraintZone::merge_zones(int num_vertices, std::vector<ConstraintZone>& zones)
{

    DisjointSets dj(num_vertices);
    for (const auto& zone : zones) {
        size_t num_stencil = zone.stencil.size();
        for (size_t i = 1; i < num_stencil; ++i) {
            dj.make_union(zone.stencil[0], zone.stencil[i]);
        }
    }

    std::vector<ConstraintZone> old_zones;
    std::swap(zones, old_zones);
    zones.reserve(old_zones.size());
    std::vector<int> parent_to_zone_index;
    parent_to_zone_index.reserve(old_zones.size());

    // Loop zones and merge
    for (auto& zone : old_zones) {

        int parent = dj.find(zone.stencil[0]);
        while (parent >= parent_to_zone_index.size()) {
            parent_to_zone_index.emplace_back(-1);
        }

        // new zone
        int zone_index = parent_to_zone_index[parent];
        if (zone_index < 0) {
            int new_zone_index = zones.size();
            parent_to_zone_index[parent] = new_zone_index;
            zones.emplace_back(zone);
            zones.back().index = new_zone_index;
        }
        // Existing zone: merge
        else {
            // Should I cache zones and merge all at once?
            zones[zone_index].merge(zone);
        }
    }
}

} // end ns mcl

#endif // MCL_GEOM_CONSTRAINT_ZONE_HPP