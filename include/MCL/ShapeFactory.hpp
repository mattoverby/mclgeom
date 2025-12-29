// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_GEOM_SHAPEFACTORY_HPP
#define MCL_GEOM_SHAPEFACTORY_HPP

#include <Eigen/Dense>
#include <array>
#include <map>
#include <vector>

namespace mcl {

/// @brief Makes a 3D box between lower and upper corners
template<typename T, typename DerivedV, typename DerivedF>
inline void
make_tri_box(const Eigen::Vector3<T>& bmin, const Eigen::Vector3<T>& bmax, DerivedV& V, DerivedF& F);

/// @brief Makes a 3D quad between lower and upper corners
template<typename T, typename DerivedV, typename DerivedF>
inline void
make_quad(const Eigen::Vector2<T>& bottom_left,
          const Eigen::Vector2<T>& upper_right,
          int tessellation,
          DerivedV& V,
          DerivedF& F);

/// @brief Makes a 3D sphere
template<typename T, typename DerivedV, typename DerivedF>
inline void
make_sphere(const Eigen::Vector3<T>& center, T radius, int subdivisions, DerivedV& V, DerivedF& F);

//
// Implementation
//

template<typename T, typename DerivedV, typename DerivedF>
inline void
make_tri_box(const Eigen::Vector3<T>& bmin, const Eigen::Vector3<T>& bmax, DerivedV& V, DerivedF& F)
{
    V.resize(8, 3);
    V << bmin[0], bmin[1], bmin[2], // 0
        bmax[0], bmin[1], bmin[2],  // 1
        bmax[0], bmax[1], bmin[2],  // 2
        bmin[0], bmax[1], bmin[2],  // 3
        bmin[0], bmin[1], bmax[2],  // 4
        bmax[0], bmin[1], bmax[2],  // 5
        bmax[0], bmax[1], bmax[2],  // 6
        bmin[0], bmax[1], bmax[2];  // 7
    F.resize(12, 3);
    F << 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 7, 3, 0, 4, 7, 1, 2, 6, 1, 6, 5, 3, 7, 6, 3, 6, 2, 0, 1, 5, 0, 5, 4;
}

template<typename T, typename DerivedV, typename DerivedF>
inline void
make_quad(const Eigen::Vector2<T>& bottom_left,
          const Eigen::Vector2<T>& upper_right,
          int tessellation,
          DerivedV& V,
          DerivedF& F)
{
    V.resize((tessellation + 1) * (tessellation + 1), 3);
    F.resize(tessellation * tessellation * 2, 3);
    Eigen::Vector2<T> size = upper_right - bottom_left;
    Eigen::Vector2<T> step = size / static_cast<T>(tessellation);

    // Generate vertices
    int vertex_count = 0;
    for (int j = 0; j <= tessellation; ++j) {
        for (int i = 0; i <= tessellation; ++i) {
            Eigen::Vector2<T> p2 = bottom_left + Eigen::Vector2<T>(i * step[0], j * step[1]);
            V.row(vertex_count++) = Eigen::RowVector3<T>(p2[0], p2[1], 0);
        }
    }

    auto idx = [tessellation](int i, int j) { return j * (tessellation + 1) + i; };

    // Generate faces (two triangles per quad)
    int face_count = 0;
    for (int j = 0; j < tessellation; ++j) {
        for (int i = 0; i < tessellation; ++i) {
            int i0 = idx(i, j);
            int i1 = idx(i + 1, j);
            int i2 = idx(i, j + 1);
            int i3 = idx(i + 1, j + 1);
            F.row(face_count++) = Eigen::RowVector3i(i0, i1, i3);
            F.row(face_count++) = Eigen::RowVector3i(i0, i3, i2);
        }
    }
}

template<typename T, typename DerivedV, typename DerivedF>
inline void
make_sphere(const Eigen::Vector3<T>& center, T radius, int subdivisions, DerivedV& V, DerivedF& F)
{
    auto get_midpoint =
        [](int i0, int i1, std::vector<Eigen::Vector3d>& vertices, std::map<std::pair<int, int>, int>& cache) -> int {
        auto key = std::minmax(i0, i1);
        auto it = cache.find(key);
        if (it != cache.end()) {
            return it->second;
        }

        Eigen::Vector3d mid = 0.5 * (vertices[i0] + vertices[i1]);
        mid.normalize();

        int index = static_cast<int>(vertices.size());
        vertices.emplace_back(mid);
        cache[key] = index;
        return index;
    };

    const double t = (1.0 + std::sqrt(5.0)) * 0.5;
    std::vector<Eigen::Vector3d> vertices = { { -1, t, 0 }, { 1, t, 0 }, { -1, -t, 0 }, { 1, -t, 0 },
                                              { 0, -1, t }, { 0, 1, t }, { 0, -1, -t }, { 0, 1, -t },
                                              { t, 0, -1 }, { t, 0, 1 }, { -t, 0, -1 }, { -t, 0, 1 } };

    for (auto& v : vertices) {
        v.normalize();
    }

    std::vector<std::array<int, 3>> triangles = { { 0, 11, 5 },  { 0, 5, 1 },  { 0, 1, 7 },  { 0, 7, 10 },
                                                  { 0, 10, 11 }, { 1, 5, 9 },  { 5, 11, 4 }, { 11, 10, 2 },
                                                  { 10, 7, 6 },  { 7, 1, 8 },  { 3, 9, 4 },  { 3, 4, 2 },
                                                  { 3, 2, 6 },   { 3, 6, 8 },  { 3, 8, 9 },  { 4, 9, 5 },
                                                  { 2, 4, 11 },  { 6, 2, 10 }, { 8, 6, 7 },  { 9, 8, 1 } };

    // Subdivide
    for (int r = 0; r < subdivisions; ++r) {
        std::map<std::pair<int, int>, int> midpoint_cache;
        std::vector<std::array<int, 3>> new_triangles;
        new_triangles.reserve(triangles.size() * 4);
        vertices.reserve(vertices.size() * 4);

        for (const auto& tri : triangles) {
            int a = get_midpoint(tri[0], tri[1], vertices, midpoint_cache);
            int b = get_midpoint(tri[1], tri[2], vertices, midpoint_cache);
            int c = get_midpoint(tri[2], tri[0], vertices, midpoint_cache);
            new_triangles.push_back({ tri[0], a, c });
            new_triangles.push_back({ tri[1], b, a });
            new_triangles.push_back({ tri[2], c, b });
            new_triangles.push_back({ a, b, c });
        }

        triangles.swap(new_triangles);
    }

    // Move vertices and scale
    V.resize(vertices.size(), 3);
    for (int i = 0; i < int(vertices.size()); ++i) {
        auto v = center + radius * vertices[i].cast<T>();
        V(i, 0) = v[0];
        V(i, 1) = v[1];
        V(i, 2) = v[2];
    }

    F.resize(triangles.size(), 3);
    for (int i = 0; i < int(triangles.size()); ++i) {
        F(i, 0) = triangles[i][0];
        F(i, 1) = triangles[i][1];
        F(i, 2) = triangles[i][2];
    }
}

} // end namespace mcl

#endif
