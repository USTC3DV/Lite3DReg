#include <filesystem>
#include <iostream>
#include <string>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <Eigen/Dense>

#include "nonrigid_spare_registration.h"
#include "rigid_fricp_registration.h"
#include "Types.h"

namespace fs = std::filesystem;
typedef OpenMesh::TriMesh_ArrayKernelT<> MeshType;
using Matrix3X = Eigen::Matrix<double, 3, Eigen::Dynamic>;

// ---------------------- 辅助函数 ----------------------
bool EnsureDirectoryExists(const std::string& path)
{
    if (!fs::exists(path)) {
        if (!fs::create_directories(path)) {
            std::cerr << "Failed to create directory: " << path << std::endl;
            return false;
        }
    }
    return true;
}

bool ReadMeshVerticesNormals(const std::string& filename, Matrix3X& vertices, Matrix3X& normals)
{
    MeshType mesh;
    OpenMesh::IO::Options opt = OpenMesh::IO::Options::VertexNormal;

    if (!OpenMesh::IO::read_mesh(mesh, filename, opt)) {
        std::cerr << "Failed to read mesh: " << filename << std::endl;
        return false;
    }

    // 如果没有法向量，自动计算
    if (!mesh.has_vertex_normals()) {
        mesh.request_face_normals();
        mesh.request_vertex_normals();
        mesh.update_normals();
    }

    size_t n_vertices = mesh.n_vertices();
    if (n_vertices == 0) {
        std::cerr << "Error: mesh has no vertices!" << std::endl;
        return false;
    }

    vertices.resize(3, n_vertices);
    normals.resize(3, n_vertices);

    for (auto v : mesh.vertices()) {
        auto point = mesh.point(v);
        auto normal = mesh.normal(v);
        vertices.col(v.idx()) << point[0], point[1], point[2];
        normals.col(v.idx()) << normal[0], normal[1], normal[2];
    }

    return true;
}

// ---------------------- main ----------------------
int main()
{
    // 定义输出路径
    std::string outpath0 = "./outpath0/";
    std::string outpath1 = "./outpath1/";
    std::string outpath2 = "./outpath2/";

    // 创建输出文件夹
    for (const auto& path : {outpath0, outpath1, outpath2}) {
        if (!EnsureDirectoryExists(path)) return -1;
    }

    Matrix3X target_points, target_normals;
    Matrix3X source_points, source_normals;
    Matrix3X deformed_points;

    // ------------------- 直接从文件注册 -------------------
    {
        RigidFricpRegistration fricp;
        NonrigidSpareRegistration spare;

        fricp.Reg("./data/target.ply","./data/source.ply", outpath0);
        spare.Reg("./data/target.obj","./data/source.obj", outpath0);
    } // fricp 和 spare 对象作用域结束，自动销毁

    // ------------------- 带参数注册 -------------------
    {
        RigidFricpRegistration fricp;
        NonrigidSpareRegistration spare;

        fricp.Paras_init(false,"",80,1e-5);
        fricp.Reg("./data/target.ply","./data/source.ply", outpath1);

        spare.Paras_init(50,1e-4,1e-5);
        spare.Reg("./data/target.obj","./data/source.obj", outpath1);
    }

    // ------------------- 使用顶点和法向量直接注册 -------------------
    if (!ReadMeshVerticesNormals("./data/target.ply", target_points, target_normals)) return -1;
    if (!ReadMeshVerticesNormals("./data/source.ply", source_points, source_normals)) return -1;

    {
        RigidFricpRegistration fricp;
        fricp.Read_data(target_points, source_points, target_normals, source_normals);
        fricp.Register();
        deformed_points = fricp.deformed_points_3X_;
        fricp.Output_data(outpath2, "fricp");
    }

    if (!ReadMeshVerticesNormals("./data/target.obj", target_points, target_normals)) return -1;
    if (!ReadMeshVerticesNormals("./data/source.obj", source_points, source_normals)) return -1;

    {
        NonrigidSpareRegistration spare;
        spare.Read_data(target_points, source_points, target_normals, source_normals);
        spare.Register();
        deformed_points = spare.deformed_points_3X_;
        spare.Output_data(outpath2, "spare");
    }

    std::cout << "Registration finished successfully!" << std::endl;
    return 0;
}
