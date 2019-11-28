#include <floattetwild/FloatTetDelaunay.h>

#include <floattetwild/Logger.hpp>

#include <igl/Timer.h>

#include <iterator>
#include <algorithm>
#include <bitset>

#include <floattetwild/Predicates.hpp>

#include <floattetwild/LocalOperations.h>
#include <floattetwild/MeshImprovement.h>
#include <floattetwild/MeshIO.hpp>



#ifdef FLOAT_TETWILD_USE_TBB
#include <tbb/task_scheduler_init.h>
#include <thread>
#endif

#include <floattetwild/Mesh.hpp>
#include <floattetwild/MeshIO.hpp>
#include <floattetwild/FloatTetDelaunay.h>
#include <floattetwild/FloatTetCutting.h>
#include <floattetwild/LocalOperations.h>
#include <floattetwild/MeshImprovement.h>
#include <floattetwild/Simplification.h>
#include <floattetwild/AABBWrapper.h>
#include <floattetwild/Statistics.h>
#include <floattetwild/TriangleInsertion.h>

#include <floattetwild/Logger.hpp>
#include <Eigen/Dense>

#include <igl/Timer.h>

#ifdef LIBIGL_WITH_TETGEN
#include <igl/copyleft/tetgen/tetrahedralize.h>
#endif

#include <geogram/basic/logger.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

#include <geogram/mesh/mesh.h>

#include<bitset>

using namespace Eigen;


namespace floatTetWild {
	namespace {
        void
        get_bb_corners(const Parameters &params, const std::vector<Vector3> &vertices, Vector3 &min, Vector3 &max) {
            min = vertices.front();
            max = vertices.front();

            for (size_t j = 0; j < vertices.size(); j++) {
                for (int i = 0; i < 3; i++) {
                    min(i) = std::min(min(i), vertices[j](i));
                    max(i) = std::max(max(i), vertices[j](i));
                }
            }

//            const Scalar dis = std::max((max - min).minCoeff() * params.box_scale, params.eps_input * 2);
            const Scalar dis = std::max(params.ideal_edge_length, params.eps_input * 2);
            for (int j = 0; j < 3; j++) {
                min[j] -= dis;
                max[j] += dis;
            }

            logger().debug("min = {} {} {}", min[0], min[1], min[2]);
            logger().debug("max = {} {} {}", max[0], max[1], max[2]);
        }

        bool comp(const std::array<int, 4> &a, const std::array<int, 4> &b) {
            return std::tuple<int, int, int>(a[0], a[1], a[2]) < std::tuple<int, int, int>(b[0], b[1], b[2]);
        }

        void match_surface_fs(Mesh &mesh, const std::vector<Vector3> &input_vertices,
                              const std::vector<Vector3i> &input_faces, std::vector<bool> &is_face_inserted) {
            std::vector<std::array<int, 4>> input_fs(input_faces.size());
            for (int i = 0; i < input_faces.size(); i++) {
                input_fs[i] = {{input_faces[i][0], input_faces[i][1], input_faces[i][2], i}};
                std::sort(input_fs[i].begin(), input_fs[i].begin() + 3);
            }
            std::sort(input_fs.begin(), input_fs.end(), comp);

//            for(auto& f: input_fs){
//                cout<<f[0]<<" "<<f[1]<<" "<<f[2]<<" "<<f[3]<<endl;
//            }
//            cout<<"/////"<<endl;

            for (auto &t: mesh.tets) {
                for (int j = 0; j < 4; j++) {
                    std::array<int, 3> f = {{t[(j + 1) % 4], t[(j + 2) % 4], t[(j + 3) % 4]}};
                    std::sort(f.begin(), f.end());
//                    cout<<f[0]<<" "<<f[1]<<" "<<f[2]<<endl;
                    auto bounds = std::equal_range(input_fs.begin(), input_fs.end(),
                                                   std::array<int, 4>({{f[0], f[1], f[2], -1}}),
                                                   comp);
//                    bool is_matched = false;
//                    int total_ori = 0;
                    for (auto it = bounds.first; it != bounds.second; ++it) {
//                        is_matched = true;
                        int f_id = (*it)[3];
                        is_face_inserted[f_id] = true;
//                        int ori = Predicates::orient_3d(mesh.tet_vertices[t[j]].pos,
//                                                        input_vertices[input_faces[f_id][0]],
//                                                        input_vertices[input_faces[f_id][1]],
//                                                        input_vertices[input_faces[f_id][2]]);
//                        if (ori == Predicates::ORI_POSITIVE)
//                            total_ori++;
//                        else if (ori == Predicates::ORI_NEGATIVE)
//                            total_ori--;
                    }
//                    if (is_matched)
//                        t.is_surface_fs[j] = total_ori;
//                    else
//                        t.is_surface_fs[j] = NOT_SURFACE;

//                    if(is_matched)
//                        cout<<"matched: "<<total_ori<<endl;
                }
            }
        }

        void match_bbox_fs(Mesh &mesh, const Vector3 &min, const Vector3 &max) {
            auto get_bbox_fs = [&](const MeshTet &t, int j) {
                std::array<int, 6> cnts = {{0, 0, 0, 0, 0, 0}};
                for (int k = 0; k < 3; k++) {
                    Vector3 &pos = mesh.tet_vertices[t[(j + k + 1) % 4]].pos;
                    for (int n = 0; n < 3; n++) {
                        if (pos[n] == min[n])
                            cnts[n * 2]++;
                        else if (pos[n] == max[n])
                            cnts[n * 2 + 1]++;
                    }
                }
                for (int i = 0; i < cnts.size(); i++) {
                    if (cnts[i] == 3)
                        return i;
                }
                return NOT_BBOX;
            };

            for (auto &t: mesh.tets) {
                for (int j = 0; j < 4; j++) {
                    t.is_bbox_fs[j] = get_bbox_fs(t, j);
                }
            }
        }

        void
        compute_voxel_points(const Vector3 &min, const Vector3 &max, const Parameters &params, const AABBWrapper &tree,
                             std::vector<Vector3> &voxels) {
            const Vector3 diag = max - min;
            Vector3i n_voxels = (diag / (params.bbox_diag_length * params.box_scale)).cast<int>();

            for (int d = 0; d < 3; ++d)
                n_voxels(d) = std::max(n_voxels(d), 1);

            const Vector3 delta = diag.array() / n_voxels.array().cast<Scalar>();

            voxels.clear();
            voxels.reserve((n_voxels(0) + 1) * (n_voxels(1) + 1) * (n_voxels(2) + 1));

//            const double sq_distg = std::min(params.ideal_edge_length / 2, 10 * params.eps);
            const double sq_distg = 10 * params.eps;
            GEO::vec3 nearest_point;

            for (int i = 0; i <= n_voxels(0); ++i) {
                const Scalar px = (i == n_voxels(0)) ? max(0) : (min(0) + delta(0) * i);
                for (int j = 0; j <= n_voxels(1); ++j) {
                    const Scalar py = (j == n_voxels(1)) ? max(1) : (min(1) + delta(1) * j);
                    for (int k = 0; k <= n_voxels(2); ++k) {
                        const Scalar pz = (k == n_voxels(2)) ? max(2) : (min(2) + delta(2) * k);

                        const GEO::vec3 gp(px, py, pz);
                        Scalar dist = sqrt(tree.project_to_sf(gp));

                        if (dist > sq_distg)
                            voxels.emplace_back(px, py, pz);
                    }
                }
            }
        }
    }

//#include <igl/unique_rows.h>
//#include <floattetwild/Predicates.hpp>
//    extern "C" floatTetWild::Scalar orient3d(const floatTetWild::Scalar *pa, const floatTetWild::Scalar *pb, const floatTetWild::Scalar *pc, const floatTetWild::Scalar *pd);

	void FloatTetDelaunay::tetrahedralize(const std::vector<Vector3>& input_vertices, const std::vector<Vector3i>& input_faces, const AABBWrapper &tree,
	        Mesh &mesh, std::vector<bool> &is_face_inserted) {
        const Parameters &params = mesh.params;
        auto &tet_vertices = mesh.tet_vertices;
        auto &tets = mesh.tets;

        is_face_inserted.resize(input_faces.size(), false);

        Vector3 min, max;
        get_bb_corners(params, input_vertices, min, max);
        mesh.params.bbox_min = min;
        mesh.params.bbox_max = max;

        std::vector<Vector3> boxpoints; //(8);
        // for (int i = 0; i < 8; i++) {
        //     auto &p = boxpoints[i];
        //     std::bitset<sizeof(int) * 8> flag(i);
        //     for (int j = 0; j < 3; j++) {
        //         if (flag.test(j))
        //             p[j] = max[j];
        //         else
        //             p[j] = min[j];
        //     }
        // }


        std::vector<Vector3> voxel_points;
        compute_voxel_points(min, max, params, tree, voxel_points);

        const int n_pts = input_vertices.size() + boxpoints.size() + voxel_points.size();
        tet_vertices.resize(n_pts);
//        std::vector<double> V_d;
//        V_d.resize(n_pts * 3);

        size_t index = 0;
        int offset = 0;
        for (int i = 0; i < input_vertices.size(); i++) {
            tet_vertices[offset + i].pos = input_vertices[i];
            // tet_vertices[offset + i].is_on_surface = true;
//            for (int j = 0; j < 3; j++)
//                V_d[index++] = input_vertices[i](j);
        }
        offset += input_vertices.size();
        for (int i = 0; i < boxpoints.size(); i++) {
            tet_vertices[i + offset].pos = boxpoints[i];
            // tet_vertices[i + offset].is_on_bbox = true;
//            for (int j = 0; j < 3; j++)
//                V_d[index++] = boxpoints[i](j);
        }
        offset += boxpoints.size();
        for (int i = 0; i < voxel_points.size(); i++) {
            tet_vertices[i + offset].pos = voxel_points[i];
            // tet_vertices[i + offset].is_on_bbox = false;
//            for (int j = 0; j < 3; j++)
//                V_d[index++] = voxel_points[i](j);
        }

        std::vector<double> V_d;
        V_d.resize(n_pts * 3);
        for (int i = 0; i < tet_vertices.size(); i++) {
            for (int j = 0; j < 3; j++)
                V_d[i * 3 + j] = tet_vertices[i].pos[j];
        }

        GEO::Delaunay::initialize();
        GEO::Delaunay_var T = GEO::Delaunay::create(3, "BDEL");
        T->set_vertices(n_pts, V_d.data());
        //
        tets.resize(T->nb_cells());
        const auto &tet2v = T->cell_to_v();
        for (int i = 0; i < T->nb_cells(); i++) {
            for (int j = 0; j < 4; ++j) {
                const int v_id = tet2v[i * 4 + j];

                tets[i][j] = v_id;
                tet_vertices[v_id].conn_tets.push_back(i);
            }
            std::swap(tets[i][1], tets[i][3]);
        }

        for (int i = 0; i < mesh.tets.size(); i++) {
            auto &t = mesh.tets[i];
            if (is_inverted(mesh.tet_vertices[t[0]].pos, mesh.tet_vertices[t[1]].pos,
                            mesh.tet_vertices[t[2]].pos, mesh.tet_vertices[t[3]].pos)) {
                cout << "EXIT_INV" << endl;
                exit(0);
            }
        }
        match_bbox_fs(mesh, min, max);
    }

    void FloatTetDelaunay::fTetWildTetrahedralize(floatTetWild::Parameters& params, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::Matrix<Scalar, Eigen::Dynamic, 3>& surface_verts, Eigen::Matrix<int, Eigen::Dynamic, 3>& surface_faces) {
        GEO::initialize();
    //    exactinit();

        std::vector<int> indices(20);
        std::iota(std::begin(indices), std::end(indices), 0);
        floatTetWild::Random::shuffle(indices);
        for (int a : indices)
            std::cout << a << " ";
        std::cout << std::endl;

        // Import standard command line arguments, and custom ones
        GEO::CmdLine::import_arg_group("standard");
        GEO::CmdLine::import_arg_group("pre");
        GEO::CmdLine::import_arg_group("algo");

        bool run_tet_gen = false;
        bool skip_simplify = false;

        Mesh mesh;
        mesh.params = params;
    int boolean_op = -1;
    unsigned int max_threads = std::numeric_limits<unsigned int>::max();
    #ifdef FLOAT_TETWILD_USE_TBB
        std::cout << "USING TBB\n";
        const size_t MB = 1024 * 1024;
        const size_t stack_size = 64 * MB;
        unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
        num_threads = std::min(max_threads, num_threads);
        params.num_threads = num_threads;
        std::cout << "TBB threads " << num_threads << std::endl;
        tbb::task_scheduler_init scheduler(num_threads, stack_size);
    #else
        std::cout << "NOT USING TBB\n";
    #endif

    //    if(params.is_quiet){
    //        std::streambuf *orig_buf = cout.rdbuf();
    //        cout.rdbuf(NULL);
    //    }

        Logger::init(!params.is_quiet, params.log_path);
        params.log_level = std::max(0, std::min(6, params.log_level));
        spdlog::set_level(static_cast<spdlog::level::level_enum>(params.log_level));
        spdlog::flush_every(std::chrono::seconds(3));

        GEO::Logger *geo_logger = GEO::Logger::instance();
        geo_logger->unregister_all_clients();
        // geo_logger->register_client(new GeoLoggerForward(logger().clone("geogram")));
        geo_logger->set_pretty(false);


        if (params.output_path.empty())
            params.output_path = params.input_path;
        if (params.log_path.empty())
            params.log_path = params.output_path;


        std::string output_mesh_name = params.output_path;
        if (params.output_path.size() > 3
            && params.output_path.substr(params.output_path.size() - 3, params.output_path.size()) == "msh")
            output_mesh_name = params.output_path;
        else if (params.output_path.size() > 4
                && params.output_path.substr(params.output_path.size() - 4, params.output_path.size()) == "mesh")
            output_mesh_name = params.output_path;
        else
            output_mesh_name = params.output_path + "_" + params.postfix + ".msh";

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        std::vector<Vector3> input_vertices;
        std::vector<Vector3i> input_faces;
        std::vector<int> input_tags;

        if (!params.tag_path.empty()) {
            input_tags.reserve(input_faces.size());
            std::string line;
            std::ifstream fin(params.tag_path);
            if (fin.is_open()) {
                while (getline(fin, line)) {
                    input_tags.push_back(std::stoi(line));
                }
                fin.close();
            }
        }


        igl::Timer timer;

        GEO::Mesh sf_mesh;
        if (!MeshIO::load_mesh(V, F, input_vertices, input_faces, sf_mesh, input_tags)) {
            logger().error("Unable to load mesh at {}", params.input_path);
        } else if (input_vertices.empty() || input_faces.empty()) {
        }

        if (input_tags.size() != input_faces.size()) {
            input_tags.resize(input_faces.size());
            std::fill(input_tags.begin(), input_tags.end(), 0);
        }
        AABBWrapper tree(sf_mesh);

        // if (!params.init(tree.get_sf_diag())) {
        // }

    #ifdef LIBIGL_WITH_TETGEN
        if(run_tet_gen)
        {
            Eigen::MatrixXd tetgen_pts(input_vertices.size(), 3);
            Eigen::MatrixXi tetgen_faces(input_faces.size(), 3);

            for(size_t i = 0; i < input_vertices.size(); ++i)
            {
                tetgen_pts.row(i) = input_vertices[i].cast<double>();
            }

            for(size_t i = 0; i < input_faces.size(); ++i)
            {
                tetgen_faces.row(i) = input_faces[i];
            }

            std::stringstream buf;
            buf.precision(100);
            buf.setf(std::ios::fixed, std::ios::floatfield);
            buf<<"Qpq2.0a"<<params.ideal_edge_length*params.ideal_edge_length*params.ideal_edge_length*sqrt(2.)/12.;

            Eigen::MatrixXi tetgen_generated_tets;
            Eigen::MatrixXd tetgen_generated_points;
            Eigen::MatrixXi tetgen_generated_faces;

            timer.start();
            igl::copyleft::tetgen::tetrahedralize(tetgen_pts, tetgen_faces, buf.str(), tetgen_generated_points, tetgen_generated_tets, tetgen_generated_faces);
            timer.stop();
            logger().info("Tetgen time {}s", timer.getElapsedTimeInSec());
            stats().record(StateInfo::tetgen_id, timer.getElapsedTimeInSec(), tetgen_generated_points.rows(), tetgen_generated_tets.rows(), 0, 0);
        }
    #endif

        stats().record(StateInfo::init_id, 0, input_vertices.size(), input_faces.size(), -1, -1);

        timer.start();
        simplify(input_vertices, input_faces, input_tags, tree, params, skip_simplify);
        tree.init_b_mesh_and_tree(input_vertices, input_faces);
        logger().info("preprocessing {}s", timer.getElapsedTimeInSec());
        logger().info("");
        stats().record(StateInfo::preprocessing_id, timer.getElapsedTimeInSec(), input_vertices.size(),
                    input_faces.size(), -1, -1);
        if (params.log_level <= 1)
            output_component(input_vertices, input_faces, input_tags);

        timer.start();
        std::vector<bool> is_face_inserted(input_faces.size(), false);
        FloatTetDelaunay::tetrahedralize(input_vertices, input_faces, tree, mesh, is_face_inserted);
        logger().info("#v = {}", mesh.get_v_num());
        logger().info("#t = {}", mesh.get_t_num());
        logger().info("tetrahedralizing {}s", timer.getElapsedTimeInSec());
        logger().info("");
        stats().record(StateInfo::tetrahedralization_id, timer.getElapsedTimeInSec(), mesh.get_v_num(), mesh.get_t_num(),
                    -1, -1);

        timer.start();
        insert_triangles(input_vertices, input_faces, input_tags, mesh, is_face_inserted, tree, false);
        logger().info("cutting {}s", timer.getElapsedTimeInSec());
        logger().info("");
        stats().record(StateInfo::cutting_id, timer.getElapsedTimeInSec(), mesh.get_v_num(), mesh.get_t_num(),
                    mesh.get_max_energy(), mesh.get_avg_energy(),
                    std::count(is_face_inserted.begin(), is_face_inserted.end(), false));

    //    timer.start();
    ////    cutting(input_vertices, input_faces, mesh, is_face_inserted, tree);
    //    cutting(input_vertices, input_faces, input_tags, mesh, is_face_inserted, tree);
    //    logger().info("cutting {}s", timer.getElapsedTimeInSec());
    //    logger().info("");
    //    stats().record(StateInfo::cutting_id, timer.getElapsedTimeInSec(), mesh.get_v_num(), mesh.get_t_num(),
    //                                                   mesh.get_max_energy(), mesh.get_avg_energy(),
    //                                                   std::count(is_face_inserted.begin(), is_face_inserted.end(), false));

        timer.start();
        optimization(input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree, {{1, 1, 1, 1}});
        logger().info("mesh optimization {}s", timer.getElapsedTimeInSec());
        logger().info("");
        stats().record(StateInfo::optimization_id, timer.getElapsedTimeInSec(), mesh.get_v_num(), mesh.get_t_num(),
                    mesh.get_max_energy(), mesh.get_avg_energy());

        timer.start();
        correct_tracked_surface_orientation(mesh, tree);
        logger().info("correct_tracked_surface_orientation done");
        if (boolean_op < 0) {
            if (params.smooth_open_boundary) {
                smooth_open_boundary(mesh, tree);
                for (auto &t: mesh.tets) {
                    if (t.is_outside)
                        t.is_removed = true;
                }
            } else
                filter_outside(mesh);
        } else
            boolean_operation(mesh, boolean_op);
        if(params.manifold_surface){
    //        MeshIO::write_mesh(params.output_path + "_" + params.postfix + "_non_manifold.msh", mesh, false);
            manifold_surface(mesh);
        }
        stats().record(StateInfo::wn_id, timer.getElapsedTimeInSec(), mesh.get_v_num(), mesh.get_t_num(),
                    mesh.get_max_energy(), mesh.get_avg_energy());
        logger().info("after winding number");
        logger().info("#v = {}", mesh.get_v_num());
        logger().info("#t = {}", mesh.get_t_num());
        logger().info("winding number {}s", timer.getElapsedTimeInSec());
        logger().info("");


    //    if (params.output_path.size() > 3
    //        && params.output_path.substr(params.output_path.size() - 3, params.output_path.size()) == "msh")
    //        MeshIO::write_mesh(params.output_path, mesh, false);
    //    else if (params.output_path.size() > 4
    //             && params.output_path.substr(params.output_path.size() - 4, params.output_path.size()) == "mesh")
    //        MeshIO::write_mesh(params.output_path, mesh, false);
    //    else
    //        MeshIO::write_mesh(params.output_path + "_" + params.postfix + ".msh", mesh, false);

        //fortest
        std::vector<Scalar> colors(mesh.tets.size(), -1);
        for (int i = 0; i < mesh.tets.size(); i++) {
            if (mesh.tets[i].is_removed)
                continue;
            colors[i] = mesh.tets[i].quality;
        }
        //fortest
        // MeshIO::write_mesh(output_mesh_name, mesh, false, colors);
        // MeshIO::write_surface_mesh(params.output_path + "_" + params.postfix + "_sf.obj", mesh, false);


        
            const auto skip_tet    = [&mesh](const int i) { return mesh.tets[i].is_removed; };
            const auto skip_vertex = [&mesh](const int i) { return mesh.tet_vertices[i].is_removed; };
            MeshIO::extract_surface_mesh(mesh, skip_tet, skip_vertex, surface_verts, surface_faces);

        // igl::write_triangle_mesh(path, V_sf, F_sf);

        timer.stop();
        logger().info(" took {}s", timer.getElapsedTime());
    }

}
