// Copyright Matt Overby 2025.
// Distributed under the MIT License.

#ifndef MCL_GEOM_VISUALDEBUG_HPP
#define MCL_GEOM_VISUALDEBUG_HPP 1

#include <Eigen/Dense>
#include <mutex>
#include <vector>

namespace mcl {

/// @brief A tool for visual debugging: cache geometry to draw on next frame.
///
/// Example:
///   // include viewer before VisualDebug.hpp to access "set_data_igl" function.
///   #include <igl/opengl/glfw/Viewer.h>
///   #include <MCL/VisualDebug.hpp>
///   ...
///   mcl::VisualDebug &vd = mcl::VisualDebug::get_instance(); // singleton
///   vd.add_point(Eigen::Vector3d(0,0,0), mcl::VisualDebug::Red); // thread-safe (mutex)
///   igl::opengl::glfw::Viewer viewer;
///   viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer&) -> bool {
///        viewer.data(0).set_mesh(V, F);
///        vd.set_data_igl(viewer, 1, false);
///        return false;
///    };
///    viewer.launch();
///
class VisualDebug
{
  public:
    static VisualDebug& get_instance()
    {
        static VisualDebug instance;
        return instance;
    }

    inline static const Eigen::Vector3d Red{ 1.0, 0.0, 0.0 };
    inline static const Eigen::Vector3d Green{ 0, 1, 0 };
    inline static const Eigen::Vector3d Blue{ 0, 0, 1 };
    inline static const Eigen::Vector3d Yellow{ 1, 1, 0 };
    inline static const Eigen::Vector3d DarkGray{ 0.4, 0.4, 0.4 };
    inline static const Eigen::Vector3d Gray{ 0.4, 0.4, 0.4 };

    VisualDebug(VisualDebug const&) = delete;

    VisualDebug operator=(VisualDebug const&) = delete;

#ifdef IGL_OPENGL_GLFW_VIEWER_H
    /// @brief Sets cache to igl mesh at specified index. Clears existings points and lines.
    /// Optionally clear buffered data.
    void set_data_igl(igl::opengl::glfw::Viewer& viewer, size_t mesh_index, bool clear_cache = true);
#endif

    /// @brief Buffer a point
    template<typename VectorType>
    void add_point(const VectorType& p, const Eigen::Vector3d& c = Red);

    /// @brief Buffer multiple points; color defaults to red
    template<typename MatrixType>
    void add_points(const MatrixType& p, const MatrixType& c = MatrixType());

    /// @brief Buffer a line segment (p0 to p1)
    template<typename VectorType>
    void add_line(const VectorType& p0, const VectorType& p1, const Eigen::Vector3d& c = Red);

    /// @brief Buffer multiple line segments (p0 to p1)
    /// color defaults to red
    template<typename MatrixType>
    void add_lines(const MatrixType& p0, const MatrixType& p1, const MatrixType& c = MatrixType());

    /// @brief Buffer a wireframe box
    template<typename VectorType>
    void add_wireframe_box(const VectorType& min, const VectorType& max, const Eigen::Vector3d& c = Red);

    /// @brief Clears all buffered data.
    void clear();

  protected:
    VisualDebug() {}

    virtual ~VisualDebug() {}

    std::mutex mtx;
    std::vector<Eigen::Vector3d> pts;
    std::vector<Eigen::Vector3d> pt_colors;
    std::vector<Eigen::Vector3d> lines0, lines1;
    std::vector<Eigen::Vector4d> line_colors;

    template<typename VectorType>
    Eigen::Vector3d to_vec3d(const VectorType& p)
    {
        return Eigen::Vector3d(p[0], p[1], (p.size() >= 3 ? double(p[2]) : 0.0));
    }

    template<typename VectorVectorType>
    Eigen::MatrixXd to_matrixXd(VectorVectorType& vec)
    {
        Eigen::MatrixXd result;
        if (vec.empty()) {
            return result;
        }
        result.resize(vec.size(), 3);
        for (size_t i = 0; i < vec.size(); ++i) {
            result.row(i) = vec[i].eval();
        }
        return result;
    }

}; // class VisualDebug

//
// Implementation
//

#ifdef IGL_OPENGL_GLFW_VIEWER_H
void
VisualDebug::set_data_igl(igl::opengl::glfw::Viewer& viewer, size_t mesh_index, bool clear_cache)
{
    std::lock_guard<std::mutex> lock(mtx);

    // Create mesh index if we don't have it.
    while (mesh_index >= viewer.data_list.size()) {
        bool visible = mesh_index == int(viewer.data_list.size());
        viewer.append_mesh(visible);
    }

    // Copy to matrix and set igl

    if (!pts.empty() && pts.size() == pt_colors.size()) {
        viewer.data(mesh_index).clear_points();
        Eigen::MatrixXd P = to_matrixXd(pts);
        Eigen::MatrixXd PC = to_matrixXd(pt_colors);
        viewer.data(mesh_index).add_points(P, PC);
    }

    if (!lines0.empty() && lines0.size() == lines1.size() && lines0.size() == line_colors.size()) {
        viewer.data(mesh_index).clear_edges();
        Eigen::MatrixXd L0 = to_matrixXd(lines0);
        Eigen::MatrixXd L1 = to_matrixXd(lines1);
        Eigen::MatrixXd LC = to_matrixXd(line_colors);
        viewer.data(mesh_index).add_edges(L0, L1, LC);
    }

    if (clear_cache) {
        clear();
    }
}
#endif

template<typename VectorType>
void
VisualDebug::add_point(const VectorType& p, const Eigen::Vector3d& c)
{
    std::lock_guard<std::mutex> lock(mtx);
    pts.emplace_back(to_vec3d(p));
    pt_colors.emplace_back(c);
}

template<typename MatrixType>
void
VisualDebug::add_points(const MatrixType& p, const MatrixType& c)
{
    std::lock_guard<std::mutex> lock(mtx);
    pts.reserve(pts.size() + p.rows());
    pt_colors.reserve(pt_colors.size() + p.rows());
    for (size_t i = 0; i < p.rows(); ++i) {
        pts.emplace_back(to_vec3d(p.row(i)));
        Eigen::Vector3i ci = i < c.rows() ? to_vec3d(c.row(i)) : Red;
        pt_colors.emplace_back(ci);
    }
}

template<typename VectorType>
void
VisualDebug::add_line(const VectorType& p0, const VectorType& p1, const Eigen::Vector3d& c)
{
    std::lock_guard<std::mutex> lock(mtx);
    lines0.emplace_back(to_vec3d(p0));
    lines1.emplace_back(to_vec3d(p1));
    line_colors.emplace_back(c);
}

template<typename MatrixType>
void
VisualDebug::add_lines(const MatrixType& p0, const MatrixType& p1, const MatrixType& c)
{
    std::lock_guard<std::mutex> lock(mtx);
    size_t nl = std::min(p0.size(), p1.size());
    lines0.reserve(lines0.size() + nl);
    lines1.reserve(lines1.size() + nl);
    line_colors.reserve(line_colors.size() + nl);
    for (size_t i = 0; i < nl; ++i) {
        lines0.emplace_back(to_vec3d(p0.row(i)));
        lines1.emplace_back(to_vec3d(p1.row(i)));
        Eigen::Vector3i ci = i < c.rows() ? to_vec3d(c.row(i)) : Red;
        line_colors.emplace_back(ci);
    }
}

template<typename VectorType>
void
VisualDebug::add_wireframe_box(const VectorType& min_, const VectorType& max_, const Eigen::Vector3d& color)
{
    Eigen::Vector3d min = to_vec3d(min_);
    Eigen::Vector3d max = to_vec3d(max_);

    // Bottom quad
    Eigen::Vector3d a = min;
    Eigen::Vector3d b(max[0], min[1], min[2]);
    Eigen::Vector3d c(max[0], min[1], max[2]);
    Eigen::Vector3d d(min[0], min[1], max[2]);

    // Top quad
    Eigen::Vector3d e(min[0], max[1], min[2]);
    Eigen::Vector3d f(max[0], max[1], min[2]);
    Eigen::Vector3d g = max;
    Eigen::Vector3d h(min[0], max[1], max[2]);

    // bottom
    lines0.emplace_back(a);
    lines1.emplace_back(b);
    line_colors.emplace_back(color);
    lines0.emplace_back(a);
    lines1.emplace_back(d);
    line_colors.emplace_back(color);
    lines0.emplace_back(c);
    lines1.emplace_back(b);
    line_colors.emplace_back(color);
    lines0.emplace_back(c);
    lines1.emplace_back(d);
    line_colors.emplace_back(color);

    // top
    lines0.emplace_back(e);
    lines1.emplace_back(f);
    line_colors.emplace_back(color);
    lines0.emplace_back(e);
    lines1.emplace_back(h);
    line_colors.emplace_back(color);
    lines0.emplace_back(g);
    lines1.emplace_back(f);
    line_colors.emplace_back(color);
    lines0.emplace_back(g);
    lines1.emplace_back(h);
    line_colors.emplace_back(color);

    // columns
    lines0.emplace_back(d);
    lines1.emplace_back(h);
    line_colors.emplace_back(color);
    lines0.emplace_back(min);
    lines1.emplace_back(e);
    line_colors.emplace_back(color);
    lines0.emplace_back(b);
    lines1.emplace_back(f);
    line_colors.emplace_back(color);
    lines0.emplace_back(c);
    lines1.emplace_back(max);
    line_colors.emplace_back(color);
}

void
VisualDebug::clear()
{
    pts.clear();
    pt_colors.clear();
    lines0.clear();
    lines1.clear();
    line_colors.clear();
}

} // end namespace mcl

#endif // MCL_GEOM_VISUALDEBUG_HPP
