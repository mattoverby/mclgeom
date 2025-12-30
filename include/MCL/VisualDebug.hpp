// Copyright Matt Overby 2025.
// Distributed under the MIT License.

#ifndef MCL_GEOM_VISUALDEBUG_HPP
#define MCL_GEOM_VISUALDEBUG_HPP 1

#include <Eigen/Dense>
#include <vector>
#include <mutex>
#include <igl/opengl/glfw/Viewer.h>

namespace mcl
{

/// @brief A tool for visual debugging: cache geometry to draw on next frame.
///
/// Example:
///   mcl::VisualDebug &vd = mcl::VisualDebug::get_instance(); // singleton
///   vd.add_point(Eigen::Vector3d(0,0,0), mcl::VisualDebug::Red); // thread-safe (mutex)
///   igl::opengl::glfw::Viewer viewer;
///   viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer&) -> bool {
///        viewer.data(0).clear();
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

	inline static const Eigen::Vector3d Red{1.0, 0.0, 0.0};
	inline static const Eigen::Vector3d Green{0,1,0};
	inline static const Eigen::Vector3d Blue{0,0,1};
	inline static const Eigen::Vector3d Yellow{1,1,0};
	inline static const Eigen::Vector3d DarkGray{0.4,0.4,0.4};
	inline static const Eigen::Vector3d Gray{0.4,0.4,0.4};

	VisualDebug(VisualDebug const &) = delete;

	VisualDebug operator=(VisualDebug const &) = delete;

	/// @brief Sets cache to igl mesh at specified index. Clears existings points and lines.
	/// Optionally clear buffered data.
	void set_data_igl(igl::opengl::glfw::Viewer &viewer, int mesh_index, bool clear_cache = true);

	/// @brief Buffer a point
	template<typename VectorType>
	void add_point(
		const VectorType &p,
		const Eigen::Vector3d &c=Red);
#if 0
	/// @brief Buffer multiple points
	template<int dim>
	void add_points(
		const Eigen::VectorXd &p,
		const Eigen::Vector4d &c=Red);

	/// @brief Buffer a line segment (p0 to p1)
	template<int dim>
	void add_line(
		const Eigen::Matrix<double,dim,1> &p0,
		const Eigen::Matrix<double,dim,1> &p1,
		const Eigen::Vector4d &c=Red);

	/// @brief Buffer multiple line segments (p0 to p1)
	template<int dim>
	void add_lines(
		const Eigen::VectorXd &p0,
		const Eigen::VectorXd &p1,
		const Eigen::Vector4d &c=Red);

	/// @brief Buffer a wireframe box
	void add_wireframe_box(
		const Eigen::Vector3d &min,
		const Eigen::Vector3d &max,
		const Eigen::Vector4d &c=Red);
#endif

	/// @brief Clears all buffered data.
	void clear();


protected:
	VisualDebug() {}

	virtual ~VisualDebug() {}

	std::mutex mtx;
	std::vector<Eigen::Vector3d> pts;
	std::vector<Eigen::Vector3d> pt_colors;
#if 0
	std::vector<Eigen::Vector3d> lines0, lines1;
	std::vector<Eigen::Vector4d> line_colors;
	std::vector<Eigen::Vector3d> faces0, faces1, faces2;
	std::vector<Eigen::Vector4d> face_colors;
	std::vector<Eigen::Vector3d> arrows0, arrows1;
	std::vector<Eigen::Vector4d> arrow_colors;
	Eigen::MatrixXd mesh_vertices;
	Eigen::MatrixXd mesh_colors;
	Eigen::MatrixXi mesh_indices;

	// Here are some shapes created via constructor.
	// Every three vertices is a face.
	std::vector<Eigen::Vector3d> cylinder;
	std::vector<Eigen::Vector3d> cone;

	// For everything else (generic shapes), append them to vertices and faces
	// Called by append_faces
	void append_shapes( Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &C );

	Eigen::Vector3d clamp( const Eigen::Vector3d &a ){
		auto c01 = []( const double &aa ){ return std::min(1.0,std::max(0.0,aa)); };
		return Eigen::Vector3d( c01(a[0]), c01(a[1]), c01(a[2]) );
	}
	Eigen::Vector4d clamp( const Eigen::Vector4d &a ){
		auto c01 = []( const double &aa ){ return std::min(1.0,std::max(0.0,aa)); };
		return Eigen::Vector4d( c01(a[0]), c01(a[1]), c01(a[2]), c01(a[3]) );
	}

#endif

	template<typename VectorType>
	Eigen::Vector3d to_vec3d(const VectorType &p){
		return Eigen::Vector3d(p[0], p[1], (p.size() >= 3 ? double(p[2]) : 0.0));
	}

	template<typename VectorVectorType>
	Eigen::MatrixXd to_matrixXd(VectorVectorType &vec)
	{
		Eigen::MatrixXd result;
		if (vec.empty()) {
			return result;
		}
		result.resize(vec.size(), 3);
		for (size_t i=0; i<vec.size(); ++i) {
			result.row(i) = vec[i].eval();
		}
		return result;
	}

}; // class VisualDebug

//
// Implementation
//

void VisualDebug::set_data_igl(igl::opengl::glfw::Viewer &viewer, int mesh_index, bool clear_cache)
{
	std::lock_guard<std::mutex> lock(mtx);

	// Create mesh index if we don't have it.
	while(mesh_index >= int(viewer.data_list.size()))
	{
		bool visible = mesh_index == int(viewer.data_list.size());
		viewer.append_mesh(visible);
	}

	viewer.data(mesh_index).clear_points();
	viewer.data(mesh_index).clear_edges();

	// Copy to matrix and set igl
	if (!pts.empty() && pts.size() == pt_colors.size()) {
		Eigen::MatrixXd P = to_matrixXd(pts);
		Eigen::MatrixXd PC = to_matrixXd(pt_colors);
		viewer.data(mesh_index).set_points(P, PC);
	}

	if (clear_cache) {
		clear();
	}
}

template<typename VectorType>
void VisualDebug::add_point(const VectorType &p, const Eigen::Vector3d &c)
{
	std::lock_guard<std::mutex> lock(mtx);
	pts.emplace_back(to_vec3d(p));
	pt_colors.emplace_back(c);
}

void VisualDebug::clear()
{
	pts.clear();
	pt_colors.clear();
}

} // end namespace mcl

#endif // MCL_GEOM_VISUALDEBUG_HPP
