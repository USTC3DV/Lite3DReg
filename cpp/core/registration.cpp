#include "registration.h"
#include "io.h"
#include "tools.h"
#include <median.h>

void Registration::FindClosestPoints(VPairs & corres)
{
    corres.resize(n_src_vertex_);

    #pragma omp parallel for
    for (int i = 0; i < n_src_vertex_; i++)
    {
        Scalar mini_dist2;
        int idx = target_tree_->closest(src_mesh_->point(src_mesh_->vertex_handle(i)).data(), mini_dist2);
        Closest c;
        c.src_idx = i;
        c.position = tar_points_.col(idx);
        c.normal = Vec2Eigen(tar_mesh_->normal(tar_mesh_->vertex_handle(idx)));
        c.min_dist2 = mini_dist2;
        c.tar_idx = idx;
        corres[i] = c;
    }
}
void Registration::FindClosestPoints(KDtree* target_tree_tem,VectorX & deformed_v,VPairs & corres)
{
	corres.resize(n_src_vertex_);
#pragma omp parallel for
	for (int i = 0; i < n_src_vertex_; i++)
	{
		Scalar mini_dist2;
		int idx = target_tree_tem->closest(deformed_v.data() + 3*i, mini_dist2);
		Closest c;
		c.src_idx = i;
		c.position = tar_points_.col(idx);
		c.min_dist2 = mini_dist2;
		c.tar_idx = idx;
		corres[i] = c;
    }

}
void Registration::FindClosestPoints(KDtree* target_tree_tem,VPairs & corres, VectorX & deformed_v, std::vector<size_t>& sample_indices)
{
    corres.resize(sample_indices.size());
#pragma omp parallel for
    for(int i = 0; i < sample_indices.size(); i++)
    {
        int sidx = sample_indices[i];
        Scalar mini_dist2;
        int tidx = target_tree_tem->closest(deformed_v.data() + 3*sidx, mini_dist2);
            Closest c;
        c.src_idx = sidx;
        c.position = tar_points_.col(tidx);
        c.min_dist2 = mini_dist2;
        c.tar_idx = tidx;
        corres[i] = c;
    }
}
double Registration::FindKnearestMed(const KDtree& kdtree,
                           const Matrix3X& X, int nk)
    {
        Eigen::VectorXd X_nearest(X.cols());
#pragma omp parallel for
        for(int i = 0; i<X.cols(); i++)
        {
            int* id = new int[nk];
            double *dist = new double[nk];
            kdtree.query(X.col(i).data(), nk, id, dist);
            Eigen::VectorXd k_dist = Eigen::Map<Eigen::VectorXd>(dist, nk);
            igl::median(k_dist.tail(nk-1), X_nearest[i]);
            delete[]id;
            delete[]dist;
        }
        double med;
        igl::median(X_nearest, med);
        return sqrt(med);
    }
void Registration::LandMarkCorres(VPairs & corres)
{
    corres.clear();
    if (pars_.landmark_src.size() != pars_.landmark_tar.size())
    {
        std::cout << "Error: landmark data wrong!!" << std::endl;
    }
    n_landmark_nodes_ = pars_.landmark_tar.size();
    for (int i = 0; i < n_landmark_nodes_; i++)
    {
        Closest c;
        c.src_idx = pars_.landmark_src[i];
        OpenMesh::VertexHandle vh = tar_mesh_->vertex_handle(pars_.landmark_tar[i]);

        if (c.src_idx > n_src_vertex_ || c.src_idx < 0)
            std::cout << "Error: source index in Landmark is out of range!" << std::endl;
        if (vh.idx() < 0)
            std::cout << "Error: target index in Landmark is out of range!" << std::endl;

        c.position = Vec2Eigen(tar_mesh_->point(vh));
        c.normal = Vec2Eigen(tar_mesh_->normal(vh));
        corres.push_back(c);
	}
    std::cout << " use landmark and landmark is ... " << pars_.landmark_src.size() << std::endl;
}


bool Registration::read_landmark(const char* filename, std::vector<int>& landmark_src, std::vector<int>& landmark_tar)
{
    std::ifstream in(filename);
    std::cout << "filename = " << filename << std::endl;
    if (!in)
    {
        std::cout << "Can't open the landmark file!!" << std::endl;
        return false;
    }
    int x, y;
    landmark_src.clear();
    landmark_tar.clear();
    while (!in.eof())
    {
        if (in >> x >> y) {
            landmark_src.push_back(x);
            landmark_tar.push_back(y);
        }
    }
    in.close();
    std::cout << "landmark_src = " << landmark_src.size() << " tar = " << landmark_tar.size() << std::endl;
    return true;
}


void Registration::InitCorrespondence(VPairs & corres)
{
    if(pars_.use_landmark)
    {
        corres.clear();
        for(size_t i = 0; i < pars_.landmark_src.size(); i++)
        {
            Closest c;
            c.src_idx = pars_.landmark_src[i];
            c.tar_idx = pars_.landmark_tar[i];
            c.position = tar_points_.col(c.tar_idx);
            c.normal = Vec2Eigen(tar_mesh_->normal(tar_mesh_->vertex_handle(c.tar_idx)));
            corres.push_back(c);
        }
    }
    else
    {
        FindClosestPoints(corres);
    }
}
void Registration::Init_data()
{
    src_mesh_ = new Mesh;
    tar_mesh_ = new Mesh;
    src_mesh_ = &src_mesh;
    tar_mesh_ = &tar_mesh;
    deformed_mesh =src_mesh;
    deformed_mesh_=&deformed_mesh;
    n_src_vertex_ = src_mesh_->n_vertices();
    n_tar_vertex_ = tar_mesh_->n_vertices();
    tar_points_.resize(3, n_tar_vertex_);
    target_normals_.resize(3, n_tar_vertex_);
    for (int i = 0; i < n_tar_vertex_; i++)
    {
        auto vh = tar_mesh_->vertex_handle(i);

    // 顶点坐标
        tar_points_(0, i) = tar_mesh_->point(vh)[0];
        tar_points_(1, i) = tar_mesh_->point(vh)[1];
        tar_points_(2, i) = tar_mesh_->point(vh)[2];  

    // 顶点法向量
        target_normals_(0, i) = tar_mesh_->normal(vh)[0];
        target_normals_(1, i) = tar_mesh_->normal(vh)[1];
        target_normals_(2, i) = tar_mesh_->normal(vh)[2];
    }

        // construct kd Tree
	target_tree_ = new KDtree(tar_points_);
	src_points_.resize(3, n_src_vertex_);
	src_normals_.resize(3, n_src_vertex_);
	corres_U0_.resize(3* n_src_vertex_);
	#pragma omp parallel for
	for (int i = 0; i < n_src_vertex_; i++)
	{
		Vec3 p = src_mesh_->point(src_mesh_->vertex_handle(i));
		src_points_(0, i) = p[0];
		src_points_(1, i) = p[1];
		src_points_(2, i) = p[2];
		Vec3 n = src_mesh_->normal(src_mesh_->vertex_handle(i));
		src_normals_(0, i) = n[0];
		src_normals_(1, i) = n[1];
		src_normals_(2, i) = n[2];
	}
	deformed_normals_ = src_normals_;
    deformed_points_ = Eigen::Map<VectorX>(src_points_.data(), 3 * n_src_vertex_);
}

void Registration::Read_data(const std::string& file_target,
                       const std::string& file_source)
{
    read_by_openmesh(file_source, src_mesh);
    read_by_openmesh(file_target, tar_mesh);
    return;
}
void Registration::Read_data(const Matrix3X &target_p,const Matrix3X &source_p,const Matrix3X &target_n,const Matrix3X &source_n) 
{
    SetMeshPoints(&src_mesh, source_p, source_n);
    SetMeshPoints(&tar_mesh, target_p, target_n);
    return;
}


void Registration::Read_data(const Mesh& tar,const Mesh& src)
{
    src_mesh=src;
    tar_mesh=tar;
    return;
}
Scalar Registration::SetMeshPoints(Mesh* mesh, const Matrix3X &point, const Matrix3X &point_n)
{
    int n_vertices = point.cols();
    mesh->request_vertex_normals(); // 确保法向量可以设置

    bool use_input_normals = (point_n.cols() == n_vertices);
    if (!use_input_normals) {
        std::cout << "[Info] 法向量数据不正确或为空，正在自动计算法向量..." << std::endl;
    }

    // 判断 mesh 是否为空
    bool mesh_is_empty = (mesh->n_vertices() == 0);
    std::vector<Mesh::VertexHandle> vhandles(n_vertices);

    if (mesh_is_empty) {
        // 空 mesh：添加顶点
        for (int i = 0; i < n_vertices; i++) {
            vhandles[i] = mesh->add_vertex(Vec3(point(0, i), point(1, i), point(2, i)));
        }
    } else {
        // 已有 mesh：直接覆盖顶点位置
        if (mesh->n_vertices() != n_vertices) {
            std::cerr << "[Warning] mesh 顶点数量和要替换的mesh不匹配" << std::endl;
        }
        int idx = 0;
        for (auto v : mesh->vertices()) {
            if (idx >= n_vertices) break;
            mesh->set_point(v, Vec3(point(0, idx), point(1, idx), point(2, idx)));
            vhandles[idx] = v;
            idx++;
        }
    }

// 设置法向量
#pragma omp parallel for
    for (int i = 0; i < n_vertices; i++) {
        if (use_input_normals) {
            Vec3 n(point_n(0, i), point_n(1, i), point_n(2, i));
            mesh->set_normal(vhandles[i], n);
        } else {
            Vec3 n(0.0, 0.0, 0.0); // 临时零向量，后续自动计算
            mesh->set_normal(vhandles[i], n);
        }
    }

    // 如果没有输入法向量，则自动计算
    if (!use_input_normals) {
        mesh->update_normals();
    }

    return 0;
}

void Registration::Output_data(const std::string& out_path,const std::string& method_name)
{
    Scalar gt_err = SetMeshPoints(deformed_mesh_, deformed_points_3X_, deformed_normals_);
    std::string out_file;
    out_file = out_path + method_name+"_res.ply";
    write_by_openmesh(out_file.c_str(), deformed_mesh, scale);
    std::cout<< "write the result to " << out_file << "\n" << std::endl;
    return;
}
