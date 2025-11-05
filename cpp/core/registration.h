#ifndef REGISTRATION_H
#define REGISTRATION_H
#include "Types.h"
#include "parameters.h"




class Registration
{
public:
    Registration() = default;
    virtual ~Registration() = default;
    

    Matrix3X                src_points_;                // (3,n)
	Matrix3X                src_normals_;               // (3,n)
    Matrix3X                src_vert_colors_;            // (3,n)
	Matrix3X                deformed_normals_;          // (3,n)
    Matrix3X                deformed_points_3X_;
	VectorX                 deformed_points_;           // (3n,)
    Matrix3X                tar_points_;               //  (3,n);
	Matrix3X                target_normals_;            // (3,n)
    Matrix3X                tar_vert_colors;           // (3,n)

    double scale=1;
    Mesh src_mesh;
    Mesh tar_mesh;
    Mesh deformed_mesh;
    Mesh* src_mesh_;
    Mesh* tar_mesh_;
    Mesh* deformed_mesh_;
    int n_src_vertex_;
    int n_tar_vertex_;
    int n_landmark_nodes_;   
    KDtree* target_tree_;                // correspondence paras
    spare::Parameters pars_;

    std::string res_trans_path_;
    MatrixXX res_trans;

    struct Closest{
        int src_idx; // vertex index from source model
        int tar_idx; // face index from target model
        Vector3 position;
        Vector3 normal;
        Scalar  min_dist2;
    };
    typedef std::vector<Closest> VPairs;
    VPairs correspondence_pairs_;

    VectorX corres_U0_;  // all correspondence points (3, n)

    void InitCorrespondence(VPairs & corres);
    void FindClosestPoints(VPairs & corres);
	void FindClosestPoints(KDtree* target_tree_tem,VectorX & deformed_v,VPairs & corres);
    void FindClosestPoints(KDtree* target_tree_tem,VPairs & corres, VectorX & deformed_v, std::vector<size_t>& sample_indices);
    double FindKnearestMed(const KDtree& kdtree,
                           const Matrix3X& X, int nk);
    void LandMarkCorres(VPairs & correspondence_pairs);
    bool read_landmark(const char* filename, std::vector<int>& landmark_src, std::vector<int>& landmark_tar);
    virtual void Read_data(const std::string& file_target,
                       const std::string& file_source) ;
    virtual void Read_data(const Matrix3X &target_p,const Matrix3X &source_p,const Matrix3X &target_n,const Matrix3X &source_n) ;
    virtual void Read_data(const Mesh& tar,const Mesh& src);

    virtual void Register(){}; 

    
    
    virtual void Init_data();
    Scalar SetMeshPoints(Mesh* mesh, const Matrix3X &point,const Matrix3X &point_n);
    virtual void Output_data(const std::string& out_path,const std::string& method_name) ;
    virtual void Reg(const std::string& file_target,
                       const std::string& file_source,
                       const std::string& out_path){};
};




#endif // REGISTRATION_H


