// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_GEOM_CONSTRAINT_ZONE_HPP
#define MCL_GEOM_CONSTRAINT_ZONE_HPP 1

#include "DisjointSets.hpp"
#include "SignedMeasure.hpp"

#include <Eigen/Core>

#include <unordered_map>
#include <unordered_set>

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
    T target_volume = 0;

    /// @brief Constructor
    VolumeConstraint() = default;

    /// @brief Constructor, scaling is computed as the surface area if not positive
    VolumeConstraint(const T* x, const Eigen::Vector<int, DIM + 1> stencil_, T scaling_ = 1);

    /// @brief Destructor
    ~VolumeConstraint() = default;

    /// @brief Returns (minimum) target volume
    T target_eval() const { return scaling * target_volume; }

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

    /// @brief Constructor
    ConstraintZone() = default;

    /// @brief Constructor
    ConstraintZone(int constraint_index, const int* sten, int stencil_size);

    /// @brief Destructor
    ~ConstraintZone() = default;

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
    // Set the scaling and target volume from the current state.
    // This should actually come from the rest state (if we have it).
    auto verts = get_verts<T, DIM, DIM + 1>(x, stencil.data());
    T abs_eval = std::abs(eval(verts));
    target_volume = T(0.1) * abs_eval;
    if (scaling < 0) {
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
ConstraintZone::merge_zones(int num_vertices, std::vector<ConstraintZone>& zones)
{
    if (zones.empty()) {
        return;
    }

    DisjointSets dj(num_vertices);
    for (const auto& zone : zones) {
        if (zone.stencil.empty()) {
            continue;
        }
        for (size_t i = 1; i < zone.stencil.size(); ++i) {
            dj.make_union(zone.stencil[0], zone.stencil[i]);
        }
    }

    std::unordered_map<int, std::vector<int>> root_to_zone_indices;
    root_to_zone_indices.reserve(zones.size());
    for (size_t i = 0; i < zones.size(); ++i) {
        const auto& zone = zones[i];
        if (zone.stencil.empty()) {
            continue;
        }
        int root = dj.find(zone.stencil[0]);
        root_to_zone_indices[root].push_back(i);
    }

    std::vector<ConstraintZone> merged_zones;
    merged_zones.reserve(root_to_zone_indices.size());
    for (auto& [root, zone_indices] : root_to_zone_indices) {
        (void)root;

        std::unordered_set<int> combined_constraints;
        std::unordered_set<int> combined_stencil;

        for (int zone_index : zone_indices) {
            const auto& zone = zones[zone_index];
            combined_constraints.insert(zone.constraints.begin(), zone.constraints.end());
            combined_stencil.insert(zone.stencil.begin(), zone.stencil.end());
        }

        ConstraintZone merged_zone;
        merged_zone.index = int(merged_zones.size());
        merged_zone.constraints.assign(combined_constraints.begin(), combined_constraints.end());
        merged_zone.stencil.assign(combined_stencil.begin(), combined_stencil.end());
        merged_zone.global_to_local.clear();
        for (size_t i = 0; i < merged_zone.stencil.size(); ++i) {
            merged_zone.global_to_local[merged_zone.stencil[i]] = int(i);
        }

        merged_zones.emplace_back(std::move(merged_zone));
    }

    zones.swap(merged_zones);
}

} // end ns mcl

#endif // MCL_GEOM_CONSTRAINT_ZONE_HPP