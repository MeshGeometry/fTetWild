#ifndef FLOATTETWILD_MESHIMPROVEMENT_H
#define FLOATTETWILD_MESHIMPROVEMENT_H

#include <floattetwild/Mesh.hpp>
#include <floattetwild/AABBWrapper.h>

namespace floatTetWild {
    void init(Mesh &mesh, AABBWrapper& tree);
    void optimization(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags, std::vector<bool> &is_face_inserted,
            Mesh &mesh, AABBWrapper& tree, const std::array<int, 4> &ops = {{1, 1, 1, 1}});
    void operation(const std::vector<Vector3> &input_vertices, const std::vector<Vector3i> &input_faces, const std::vector<int> &input_tags, std::vector<bool> &is_face_inserted,
            Mesh &mesh, AABBWrapper& tree, const std::array<int, 5> &ops = {{1, 1, 1, 1, 1}});
    bool update_scaling_field(Mesh &mesh, Scalar max_energy);

    int get_max_p(const Mesh &mesh);

    void get_tracked_surface(const Mesh& mesh, Eigen::Matrix<Scalar, Eigen::Dynamic, 3> &V, Eigen::Matrix<int, Eigen::Dynamic, 3> &F, int c_id = 0);
    void boolean_operation(Mesh& mesh, int op);
    void filter_outside(Mesh& mesh, bool invert_faces = false);
    void mark_outside(Mesh& mesh, bool invert_faces = false);

    void output_info(Mesh& mesh, const AABBWrapper& tree);
    void check_envelope(Mesh& mesh, const AABBWrapper& tree);
    void output_surface(Mesh& mesh, const std::string& filename);
}

#endif //FLOATTETWILD_MESHIMPROVEMENT_H
