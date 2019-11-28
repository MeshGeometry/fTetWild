#pragma once

#include <floattetwild/Mesh.hpp>
#include <floattetwild/Types.hpp>

namespace floatTetWild
{
	class MeshIO
	{
	public:
		static bool load_mesh(const std::string &path, std::vector<Vector3> &points, std::vector<Vector3i> &faces, GEO::Mesh& input, std::vector<int> &flags);
		static bool load_mesh(const Eigen::MatrixXd&     V,
                       const Eigen::MatrixXi&     F,
                       std::vector<Vector3>&  points,
                       std::vector<Vector3i>& faces,
                       GEO::Mesh&             input,
                       std::vector<int>&      flags);
		static void write_mesh(const std::string &path, const Mesh &mesh,
		        const bool do_filter = true, const std::vector<Scalar> &color = std::vector<Scalar>(), const bool binary = true);
		static void write_mesh(const std::string &path, const Mesh &mesh, const std::vector<int> &t_ids,
		        const bool do_filter = true, const bool binary = true);
		static void write_surface_mesh(const std::string &path, const Mesh &mesh, const bool only_interior=true);

		static void extract_volume_mesh(const Mesh &mesh, MatrixXs &V, Eigen::MatrixXi &T, bool only_interior = true);

		static void extract_surface_mesh(const Mesh&                               mesh,
                          const std::function<bool(int)>&           skip_tet,
                          const std::function<bool(int)>&           skip_vertex,
                          Eigen::Matrix<Scalar, Eigen::Dynamic, 3>& VS,
                          Eigen::Matrix<int, Eigen::Dynamic, 3>&    FS);
	};
}
