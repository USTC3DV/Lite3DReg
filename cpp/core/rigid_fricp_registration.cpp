#include "rigid_fricp_registration.h"
#include "tools.h"
#include "robust_norm.h"
#include <median.h>
#include "io.h"

    AffineMatrix3 RigidFricpRegistration::LogMatrix(const AffineMatrix3& T)
     {
        Eigen::RealSchur<AffineMatrix3> schur(T);
        AffineMatrix3 U = schur.matrixU();
        AffineMatrix3 R = schur.matrixT();
        std::vector<bool> selected(3, true);
        Matrix33 mat_B = Matrix33::Zero(3, 3);
        Matrix33 mat_V = Matrix33::Identity(3, 3);

        for (int i = 0; i < 3; i++)
        {
            if (selected[i] && fabs(R(i, i) - 1)> SAME_THRESHOLD)
            {
                int pair_second = -1;
                for (int j = i + 1; j <3; j++)
                {
                    if (fabs(R(j, j) - R(i, i)) < SAME_THRESHOLD)
                    {
                        pair_second = j;
                        selected[j] = false;
                        break;
                    }
                }
                if (pair_second > 0)
                {
                    selected[i] = false;
                    R(i, i) = R(i, i) < -1 ? -1 : R(i, i);
                    double theta = acos(R(i, i));
                    if (R(i, pair_second) < 0)
                    {
                        theta = -theta;
                    }
                    mat_B(i, pair_second) += theta;
                    mat_B(pair_second, i) += -theta;
                    mat_V(i, pair_second) += -theta / 2;
                    mat_V(pair_second, i) += theta / 2;
                    double coeff = 1 - (theta * R(i, pair_second)) / (2 * (1 - R(i, i)));
                    mat_V(i, i) += -coeff;
                    mat_V(pair_second, pair_second) += -coeff;
                }
            }
        }

        AffineMatrix3 LogTrim = AffineMatrix3::Zero();
        LogTrim.block(0, 0, 3, 3) = mat_B;
        LogTrim.block(0, 3, 3, 1) = mat_V * R.block(0, 3, 3, 1);
        AffineMatrix3 res = U * LogTrim * U.transpose();
        return res;
    }

    Affine3d RigidFricpRegistration::point_to_point(Matrix3X& X,
                         Matrix3X& Y,
                        const VectorX& w) {
        int dim = X.rows();
        /// Normalize weight vector
        Eigen::VectorXd w_normalized = w / w.sum();
        /// De-mean
        Eigen::VectorXd X_mean(dim), Y_mean(dim);
        for (int i = 0; i<dim; ++i) {
            X_mean(i) = (X.row(i).array()*w_normalized.transpose().array()).sum();
            Y_mean(i) = (Y.row(i).array()*w_normalized.transpose().array()).sum();
        }
        X.colwise() -= X_mean;
        Y.colwise() -= Y_mean;
        /// Compute transformation
        Affine3d transformation;
        MatrixXX sigma = X * w_normalized.asDiagonal() * Y.transpose();
        Eigen::JacobiSVD<MatrixXX> svd(sigma, Eigen::ComputeFullU | Eigen::ComputeFullV);
        if (svd.matrixU().determinant()*svd.matrixV().determinant() < 0.0) {
            Vector3 S = Vector3::Ones(dim); S(dim-1) = -1.0;
            transformation.linear() = svd.matrixV()*S.asDiagonal()*svd.matrixU().transpose();
        }
        else {
            transformation.linear() = svd.matrixV()*svd.matrixU().transpose();
        }
        transformation.translation() = Y_mean - transformation.linear()*X_mean;
        /// Re-apply mean
        X.colwise() += X_mean;
        Y.colwise() += Y_mean;
        /// Return transformation
        return transformation;
    }

     Eigen::Affine3d RigidFricpRegistration::point_to_plane(Eigen::Matrix3Xd& X,
                                   Eigen::Matrix3Xd& Y,
                                   const Eigen::Matrix3Xd& Norm,
                                   const Eigen::VectorXd& w,
                                   const Eigen::VectorXd& u) {
        /// Normalize weight vector
        Eigen::VectorXd w_normalized = w / w.sum();
        /// De-mean
        Eigen::Vector3d X_mean;
        for (int i = 0; i<3; ++i)
            X_mean(i) = (X.row(i).array()*w_normalized.transpose().array()).sum();
        X.colwise() -= X_mean;
        Y.colwise() -= X_mean;
        /// Prepare LHS and RHS
        Matrix66 LHS = Matrix66::Zero();
        Vector6 RHS = Vector6::Zero();
        Block33 TL = LHS.topLeftCorner<3, 3>();
        Block33 TR = LHS.topRightCorner<3, 3>();
        Block33 BR = LHS.bottomRightCorner<3, 3>();
        Eigen::MatrixXd C = Eigen::MatrixXd::Zero(3, X.cols());

#pragma omp parallel
        {
#pragma omp for
            for (int i = 0; i<X.cols(); i++) {
                C.col(i) = X.col(i).cross(Norm.col(i));
            }
#pragma omp sections nowait
            {
#pragma omp section
                for (int i = 0; i<X.cols(); i++) TL.selfadjointView<Eigen::Upper>().rankUpdate(C.col(i), w(i));
#pragma omp section
                for (int i = 0; i<X.cols(); i++) TR += (C.col(i)*Norm.col(i).transpose())*w(i);
#pragma omp section
                for (int i = 0; i<X.cols(); i++) BR.selfadjointView<Eigen::Upper>().rankUpdate(Norm.col(i), w(i));
#pragma omp section
                for (int i = 0; i<C.cols(); i++) {
                    double dist_to_plane = -((X.col(i) - Y.col(i)).dot(Norm.col(i)) - u(i))*w(i);
                    RHS.head<3>() += C.col(i)*dist_to_plane;
                    RHS.tail<3>() += Norm.col(i)*dist_to_plane;
                }
            }
        }
        LHS = LHS.selfadjointView<Eigen::Upper>();
        /// Compute transformation
        Eigen::Affine3d transformation;
        Eigen::LDLT<Matrix66> ldlt(LHS);
        RHS = ldlt.solve(RHS);
        transformation = Eigen::AngleAxisd(RHS(0), Eigen::Vector3d::UnitX()) *
                Eigen::AngleAxisd(RHS(1), Eigen::Vector3d::UnitY()) *
                Eigen::AngleAxisd(RHS(2), Eigen::Vector3d::UnitZ());
        transformation.translation() = RHS.tail<3>();

        /// Apply transformation
        /// Re-apply mean
        X.colwise() += X_mean;
        Y.colwise() += X_mean;
        transformation.translation() += X_mean - transformation.linear()*X_mean;
        /// Return transformation
        return transformation;
    }

void RigidFricpRegistration::point_to_point(Matrix3X& X, Matrix3X& Y,Matrix3X& Z, Vector3& source_mean_,
                        Vector3& target_mean_, ICP::Parameters& par){
        /// Build kd-tree
        KDtree kdtree(Y);
        /// Buffers
        n_src_vertex_ = X.cols();
        Matrix3X Q = Matrix3X::Zero(3, X.cols());
        VectorX W = VectorX::Zero(X.cols());
        deformed_points_ = VectorX::Zero(3 * n_src_vertex_);
        Affine3d T;
        if (par.use_init) 
        {
            T.matrix() = par.init_trans;
        }
        else 
        {
            T = Affine3d::Identity();
        }
        MatrixXX To1 = T.matrix();
        MatrixXX To2 = T.matrix();
        //Anderson Acc para
        AndersonAcceleration accelerator_;
        Affine3d SVD_T = T;
        double energy = .0, last_energy = std::numeric_limits<double>::max();
        //ground truth point clouds
        Matrix3X X_gt = X;
        if(par.has_groundtruth)
        {
            Vector3 temp_trans = par.gt_trans.col(3).head(3);
            X_gt.colwise() += source_mean_;
            X_gt = par.gt_trans.block(0, 0, 3, 3) * X_gt;
            X_gt.colwise() += temp_trans - target_mean_;
        }
        //output para
        std::string file_out = par.out_path;
        double gt_mse = 0.0;
        // dynamic welsch paras
        double nu1 = 1, nu2 = 1;
        Matrix3X X_deformed = T * X;
        deformed_points_ = Eigen::Map<const VectorX>(X_deformed.data(), 3 * n_src_vertex_);
        FindClosestPoints(target_tree_,deformed_points_,correspondence_pairs_);
        #pragma omp parallel for
        for (int i = 0; i<n_src_vertex_; ++i) 
        {
            W[i] =  correspondence_pairs_[i].min_dist2;
            Q.col(i) = correspondence_pairs_[i].position;
        }
        //dynamic welsch, calc k-nearest points with itself;
        nu2 = par.nu_end_k * FindKnearestMed(kdtree, Y, 7);
        double med1;
        Eigen::VectorXd W_sqrt = W.array().sqrt();
        igl::median(W_sqrt, med1);

        nu1 = par.nu_begin_k * med1;
        nu1 = nu1>nu2? nu1:nu2;
    
        //AA init
        accelerator_.init(par.anderson_m, (3 + 1) * (3 + 1), LogMatrix(T.matrix()).data());

        bool stop1 = false;
        while(!stop1)
        {
            /// run ICP
            int icp = 0;
            for (; icp<par.max_icp; ++icp)
            {
                bool accept_aa = false;
                energy = get_energy(W, nu1);
                
                    if (energy < last_energy) 
                    {
                        last_energy = energy;
                        accept_aa = true;
                    }
                    else
                    {
                        accelerator_.replace(LogMatrix(SVD_T.matrix()).data());
                        X_deformed = SVD_T * X;
                        deformed_points_ = Eigen::Map<const VectorX>(X_deformed.data(), 3 * n_src_vertex_);
                        FindClosestPoints(target_tree_,deformed_points_,correspondence_pairs_);
                    #pragma omp parallel for
                        for (int i = 0; i<n_src_vertex_; ++i) 
                            {
                                W[i] =  correspondence_pairs_[i].min_dist2;
                                Q.col(i) = correspondence_pairs_[i].position;
                            }
                        last_energy = get_energy(W, nu1);
                    }
                
                

                if(par.has_groundtruth)
                {
                    gt_mse = (T*X - X_gt).squaredNorm()/n_src_vertex_;
                }

                robust_weight(W, nu1);
                // Rotation and translation update
                Eigen::VectorXd W_sqrt = W.array().sqrt();
                T = point_to_point(X, Q, W_sqrt);

                //Anderson Acc
                SVD_T = T;
                
                AffineMatrix3 Trans = (Eigen::Map<const AffineMatrix3>(accelerator_.compute(LogMatrix(T.matrix()).data()).data(), 3+1, 3+1)).exp();
                T.linear() = Trans.block(0,0,3,3);
                T.translation() = Trans.block(0,3,3,1);
                X_deformed = T * X;
                deformed_points_ = Eigen::Map<const VectorX>(X_deformed.data(), 3 * n_src_vertex_);
                FindClosestPoints(target_tree_,deformed_points_,correspondence_pairs_);
            #pragma omp parallel for
                for (int i = 0; i<n_src_vertex_; ++i) 
                    {
                        W[i] =  correspondence_pairs_[i].min_dist2;
                        Q.col(i) = correspondence_pairs_[i].position;
                    }
                /// Stopping criteria
                double stop2 = (T.matrix() - To2).norm();
                To2 = T.matrix();
                if(stop2 < par.stop)
                {
                    break;
                }
            }
            
            stop1 = fabs(nu1 - nu2)<SAME_THRESHOLD? true: false;
            nu1 = nu1*par.nu_alpha > nu2? nu1*par.nu_alpha : nu2;
                
            accelerator_.reset(LogMatrix(T.matrix()).data());
            last_energy = std::numeric_limits<double>::max();
                
        }

        ///calc convergence energy
   
        last_energy = get_energy(W, nu1);
        Z = T * X;
        gt_mse = (Z-X_gt).squaredNorm()/n_src_vertex_;
        T.translation() += - T.rotation() * source_mean_ + target_mean_;
        Z.colwise() += target_mean_;

        ///save convergence result
        par.convergence_energy = last_energy;
        par.convergence_gt_mse = gt_mse;
        par.res_trans = T.matrix();
    }

    void RigidFricpRegistration::point_to_plane(Eigen::Matrix3Xd& X,
                        Eigen::Matrix3Xd& Y, Eigen::Matrix3Xd& norm_x, Eigen::Matrix3Xd& norm_y,
                        Eigen::Vector3d& source_mean_, Eigen::Vector3d& target_mean_,
                        ICP::Parameters &par) {
        /// Build kd-tree
        KDtree kdtree(Y);
        /// Buffers
        Eigen::Matrix3Xd Qp = Eigen::Matrix3Xd::Zero(3, X.cols());
        Eigen::Matrix3Xd Qn = Eigen::Matrix3Xd::Zero(3, X.cols());
        Eigen::VectorXd W = Eigen::VectorXd::Zero(X.cols());
        Eigen::Matrix3Xd ori_X = X;
        Affine3d T;
        if (par.use_init) T.matrix() = par.init_trans;
        else T = Affine3d::Identity();
        AffineMatrix3 To1 = T.matrix();
        X = T*X;

        Eigen::Matrix3Xd X_gt = X;
        if(par.has_groundtruth)
        {
            Eigen::Vector3d temp_trans = par.gt_trans.block(0, 3, 3, 1);
            X_gt = ori_X;
            X_gt.colwise() += source_mean_;
            X_gt = par.gt_trans.block(0, 0, 3, 3) * X_gt;
            X_gt.colwise() += temp_trans - target_mean_;
        }
      
        double gt_mse = 0.0;

        //Anderson Acc para
        AndersonAcceleration accelerator_;
        Affine3d LG_T = T;
        double energy = 0.0, prev_res = std::numeric_limits<double>::max(), res = 0.0;

        // Find closest point
#pragma omp parallel for
        for (int i = 0; i<X.cols(); ++i) {
            int id = kdtree.closest(X.col(i).data());
            Qp.col(i) = Y.col(id);
            Qn.col(i) = norm_y.col(id);
            W[i] = std::abs(Qn.col(i).transpose() * (X.col(i) - Qp.col(i)));
        }

        bool stop1 = false;
        while(!stop1)
        {
            /// ICP
            for(int icp=0; icp<par.max_icp; ++icp) {
                

                bool accept_aa = false;
                W = W.array().square();//如果后续合并逻辑，就要去掉
                energy = get_energy( W, par.p);
        
                Eigen::VectorXd test_w = (X-Qp).colwise().norm();
                if(par.has_groundtruth)
                {
                    gt_mse = (X - X_gt).squaredNorm()/X.cols();
                }

                /// Compute weights
                W = W.array().square();
                robust_weight( W, par.p);
                /// Rotation and translation update
                T = point_to_plane(X, Qp, Qn, W, Eigen::VectorXd::Zero(X.cols()))*T;
                /// Find closest point
#pragma omp parallel for
                for(int i=0; i<X.cols(); i++) {
                    X.col(i) = T * ori_X.col(i);
                    int id = kdtree.closest(X.col(i).data());
                    Qp.col(i) = Y.col(id);
                    Qn.col(i) = norm_y.col(id);
                    W[i] = std::abs(Qn.col(i).transpose() * (X.col(i) - Qp.col(i)));
                }

                /// Stopping criteria
                double stop2 = (T.matrix() - To1).norm();
                To1 = T.matrix();
                if(stop2 < par.stop) break;
            }
            stop1 = true;
        }

        par.res_trans = T.matrix();

        ///calc convergence energy
        W = (Qn.array()*(X - Qp).array()).colwise().sum().abs().transpose();
        W = W.array().square();
        energy = get_energy(W, par.p);
        gt_mse = (X - X_gt).squaredNorm() / X.cols();
        T.translation().noalias() += -T.rotation()*source_mean_ + target_mean_;
        X.colwise() += target_mean_;
        norm_x = T.rotation()*norm_x;

        ///save convergence result
        par.convergence_energy = energy;
        par.convergence_gt_mse = gt_mse;
        par.res_trans = T.matrix();
    
    }
    void RigidFricpRegistration::Paras_init(bool useinit,std::string fileinit,int maxiter,double stop)
    {
        pars.use_init = useinit;//not use the initial transformation
        file_init_ = fileinit;
        pars.max_icp=maxiter;
        pars.stop=stop;

    }
void RigidFricpRegistration::use_init_transform(bool a)
{
    if(a)
    {
        MatrixXX init_trans;
        read_transMat(init_trans, file_init_);
        init_trans.block(0, 3, 3, 1) /= scale;
        init_trans.block(0,3,3,1) += init_trans.block(0,0,3,3)*source_mean_ - target_mean_;
        pars.use_init = true;
        pars.init_trans = init_trans;
        //spars.init_trans = init_trans;
    }


}

void RigidFricpRegistration::Register()
{   
    // normalization
    scale = mesh_scaling(src_mesh, tar_mesh);
    Init_data();
    //scale=pointwise_normalize(tar_points_,src_points_,source_mean_,target_mean_);
    use_init_transform(pars.use_init);
    //--- Execute registration
    std::cout << "begin registration..." << std::endl;
    point_to_point(src_points_, tar_points_, deformed_points_3X_,source_mean_, target_mean_, pars);
    res_trans = pars.res_trans;
	std::cout << "Registration done!" <<std::endl;
    return ;
}
void RigidFricpRegistration::Reg(const std::string& file_target,
                       const std::string& file_source,
                       const std::string& out_path)
{
    Read_data(file_target, file_source);
    Register();
    Output_data(out_path,"FRICP");
    return;
}



