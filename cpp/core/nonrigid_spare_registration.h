#ifndef NONRIGIDREG_SPARE_H
#define NONRIGIDREG_SPARE_H

#include "nodeSampler.h"
#include "registration.h"




#ifdef USE_PARDISO
#include <Eigen/PardisoSupport>
#endif


class NonrigidSpareRegistration : public Registration {
public:
    NonrigidSpareRegistration();
    ~NonrigidSpareRegistration();

    // Adjusted parameters
    spare::Parameters pars_;

    // Non-rigid energy function parameters
    VectorX weight_d_;                        // robust weight for alignment α_i
    RowMajorSparseMatrix mat_A0_;             // symmetric coefficient matrix for linear equations
    VectorX vec_b_;                           // right-hand side for linear equations

    // Welsch parameters
    Scalar nu;

#ifdef USE_PARDISO
    Eigen::PardisoLDLT<RowMajorSparseMatrix, Eigen::Lower> solver_;
#else
    Eigen::SimplicialLDLT<RowMajorSparseMatrix, Eigen::Lower> solver_;
#endif

    // Sampling parameters
    int num_sample_nodes;                     // (r,) number of sample nodes
    int num_graph_edges;
    int num_edges;

    // Node storage structure
    svr::nodeSampler src_sample_nodes;

    // Variables
    VectorX X_;                               // (12r,) transformations of sample nodes

    // Alignment term
    VectorX nodes_P_;                         // (3n,) all sample nodes' coordinates
    RowMajorSparseMatrix align_coeff_PV0_;    // (3n, 12r) coefficient matrix F

    // Smooth term
    RowMajorSparseMatrix reg_coeff_B_;        // (6|E_G|, 12r) smoothness between nodes
    VectorX reg_right_D_;                     // (6|E_G|,) coordinate differences between xi and xj
    VectorX reg_cwise_weights_;               // (6|E_G|,) smooth weights

    // landmark term 
    RowMajorSparseMatrix    landmark_mul_;             // (12r, 12r) 
    VectorX                 tar_landmarks_;             // (3, l)    
    VectorX                 landmark_right_;            // (12r,)
    RowMajorSparseMatrix    landmark_mul_fine_;         // (3n, 3n) 
    VectorX                 landmark_right_fine_;       // (3n,)
    // Rotation matrix term
    VectorX nodes_R_;                         // (9r,) proj(A)
    RowMajorSparseMatrix rigid_coeff_L_;      // (12r, 12r) H
    RowMajorSparseMatrix rigid_coeff_J_;      // (9r, 12r) Y
    VectorX diff_UP_;                         // auxiliary matrix

    // ARAP term (coarse)
    VectorX arap_laplace_weights_;            // (6E,)
    Matrix3X local_rotations_;                // (3n, 3) R
    RowMajorSparseMatrix arap_coeff_;         // (6E, 12r) B
    RowMajorSparseMatrix arap_coeff_mul_;     // (12r, 12r) BᵀB
    VectorX arap_right_;                      // (6E,) L

    // ARAP term (fine)
    RowMajorSparseMatrix arap_coeff_fine_;    // (6E, 3n) B
    RowMajorSparseMatrix arap_coeff_mul_fine_;// (3n, 3n) BᵀB
    VectorX arap_right_fine_;                 // (6E,) Y

    // Point clouds
    RowMajorSparseMatrix normals_sum_;        // (n, 3n) N for alignment term

    // Sampling points & vertices relation
    std::vector<size_t> sampling_indices_;
    std::vector<int> vertex_sample_indices_;
    std::vector<int>        vertex_landmark_indices_;
    // KNN-neighbor indices for source points (if no faces)
    Eigen::MatrixXi src_knn_indices_;

    int align_sampling_num_ = 3000;

    // Weights of terms during optimization
    Scalar w_align;
    Scalar w_smo;
    Scalar optimize_w_align;
    Scalar optimize_w_smo;
    Scalar optimize_w_rot;
    Scalar optimize_w_arap;
    Scalar optimize_w_landmark;
    // Initialization & main steps
    void Initialize();
    void InitFromInput(Mesh& src_mesh, Mesh& tar_mesh, spare::Parameters& paras);
    virtual Scalar DoNonRigid();

    // Welsch weighting
    //void welsch_weight(VectorX& r, Scalar p);
    void InitWelschParam();
    void InitwithLandmark();

    // Rotation calculations
    void CalcNodeRotations();
    void InitRotations();
    void CalcLocalRotations(bool isCoarseAlign);

    // ARAP calculations
    void FullInARAPCoeff();
    void CalcARAPRight();
    void CalcARAPRightFine();
    // Normal calculations
    void CalcDeformedNormals();
    void InitNormalsSum();
    void CalcNormalsSum();
    void CalcLandmarkCoeff();
    // Optimization steps
    void PointwiseFineReg(Scalar nu1);
    void GraphCoarseReg(Scalar nu1);
void Paras_init(int iters = 30 ,double stopcoarse=1e-3,double stopfine=1e-4,bool uselandmark = false,std::vector<int> src=std::vector<int>(),std::vector<int> tar=std::vector<int>());
    
    void Register()override;
    void Reg(const std::string& file_target,
                       const std::string& file_source,
                       const std::string& out_path)override;
    std::string landmark_file;
    spare::Parameters paras;
    bool normalize = true; 

};



#endif
