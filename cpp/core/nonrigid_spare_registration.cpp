
#include "nonrigid_spare_registration.h"
#include "tools.h"
#include "robust_norm.h"
#include <median.h>

#if (__cplusplus >= 201402L) || (defined(_MSC_VER) && _MSC_VER >= 1800)
#define MAKE_UNIQUE std::make_unique
#else
#define MAKE_UNIQUE company::make_unique
#endif

/// @param Source (one 3D point per column)
/// @param Target (one 3D point per column)
/// @param Confidence weights


NonrigidSpareRegistration::NonrigidSpareRegistration() {
};

NonrigidSpareRegistration::~NonrigidSpareRegistration()
{
}

void NonrigidSpareRegistration::Initialize()
{

    //welsch weight paramete
    InitWelschParam();


    int knn_num_neighbor = 6;

    if(!pars_.use_geodesic_dist)
    {
        src_knn_indices_.resize(knn_num_neighbor, n_src_vertex_);
        KDtree* src_tree = new KDtree(src_points_);
#pragma omp parallel for
        for(int i = 0; i < n_src_vertex_; i++)
        {
            int* out_indices = new int[knn_num_neighbor+1];
            Scalar *out_dists = new Scalar[knn_num_neighbor+1];
            src_tree->query(src_points_.col(i).data(), knn_num_neighbor+1, out_indices, out_dists);
            for(int j = 0; j < knn_num_neighbor; j++)
            {
                src_knn_indices_(j, i) = out_indices[j+1];
            }
            delete[] out_indices;
            delete[] out_dists;
        }
        delete src_tree;
    }


    Scalar sample_radius;

    if(pars_.use_geodesic_dist)
        sample_radius = src_sample_nodes.SampleAndConstuct(*src_mesh_, pars_.uni_sample_radio,  src_points_);
    else
        sample_radius = src_sample_nodes.SampleAndConstuctFPS(*src_mesh_, pars_.uni_sample_radio, src_points_, src_knn_indices_, 4, 8);


    num_sample_nodes = src_sample_nodes.nodeSize();
    pars_.num_sample_nodes = num_sample_nodes;

    X_.resize(12 * num_sample_nodes); X_.setZero();
    align_coeff_PV0_.resize(3 * n_src_vertex_, 12 * num_sample_nodes);
    nodes_P_.resize(n_src_vertex_ * 3);

    nodes_R_.resize(9 * num_sample_nodes); nodes_R_.setZero();
    rigid_coeff_L_.resize(12 * num_sample_nodes, 12 * num_sample_nodes);
    rigid_coeff_J_.resize(12 * num_sample_nodes, 9 * num_sample_nodes);

    std::vector<Triplet> coeffv(4 * num_sample_nodes);
    std::vector<Triplet> coeffL(9 * num_sample_nodes);
    std::vector<Triplet> coeffJ(9 * num_sample_nodes);
    for (int i = 0; i < num_sample_nodes; i++)
    {
        // X_
        X_[12 * i] = 1.0;
        X_[12 * i + 4] = 1.0;
        X_[12 * i + 8] = 1.0;

        // nodes_R_
        nodes_R_[9 * i] = 1.0;
        nodes_R_[9 * i + 4] = 1.0;
        nodes_R_[9 * i + 8] = 1.0;

        for (int j = 0; j < 9; j++)
        {
            // rigid_coeff_L_
            coeffL[9 * i + j] = Triplet(12 * i + j, 12 * i + j, 1.0);
            // rigid_coeff_J_
            coeffJ[9 * i + j] = Triplet(12 * i + j, 9 * i + j, 1.0);
        }
    }
    rigid_coeff_L_.setFromTriplets(coeffL.begin(), coeffL.end());
    rigid_coeff_J_.setFromTriplets(coeffJ.begin(), coeffJ.end());


    // update coefficient matrices
    src_sample_nodes.initWeight(align_coeff_PV0_, nodes_P_,
                                reg_coeff_B_, reg_right_D_, reg_cwise_weights_);

    num_graph_edges = reg_cwise_weights_.rows();



    // update ARAP coeffs
    FullInARAPCoeff();

    local_rotations_.resize(3, n_src_vertex_ * 3);

    if(pars_.use_geodesic_dist)
    {
        num_edges = src_mesh_->n_halfedges();
        arap_right_.resize(3*src_mesh_->n_halfedges());
        arap_right_fine_.resize(3 * src_mesh_->n_halfedges());
    }
    else{
        num_edges = knn_num_neighbor*n_src_vertex_;
        arap_right_.resize(3*knn_num_neighbor*n_src_vertex_);
        arap_right_fine_.resize(3 * knn_num_neighbor*n_src_vertex_);
    }



    InitRotations();



    sampling_indices_.clear();

    // start points
    size_t startIndex = 0;
    sampling_indices_.push_back(startIndex);


    // FPS to get sampling points in align term

    VectorX minDistances(n_src_vertex_);
    minDistances.setConstant(std::numeric_limits<Scalar>::max());
    minDistances[startIndex] = 0;

    vertex_sample_indices_.resize(n_src_vertex_, -1);
    vertex_sample_indices_[startIndex] = 0;

    // repeat select farthest points
    while (sampling_indices_.size() < align_sampling_num_) {
// calculate the distance between each point with the sampling points set.
#pragma omp parallel for
        for (size_t i = 0; i < n_src_vertex_; ++i) {
            if(i==startIndex)
                continue;

            Scalar dist = (src_points_.col(startIndex) - src_points_.col(i)).norm();
            if(dist < minDistances[i])
                minDistances[i] = dist;
        }

        // choose farthest point
        int maxDistanceIndex;
        minDistances.maxCoeff(&maxDistanceIndex);
        minDistances[maxDistanceIndex] = 0;

        // add the farthest point into the sampling points set.
        sampling_indices_.push_back(maxDistanceIndex);
        startIndex= maxDistanceIndex;
        vertex_sample_indices_[startIndex] = sampling_indices_.size()-1;
    }
	if(pars_.use_landmark && pars_.landmark_src.size() > 0)
    {
       CalcLandmarkCoeff();
	   vertex_landmark_indices_.resize(n_src_vertex_, -1);
	   for(int i = 0; i < pars_.landmark_src.size(); i++)
	   {
			vertex_landmark_indices_[pars_.landmark_src[i]] = i; 
	   }
   	}
}
void NonrigidSpareRegistration::CalcLandmarkCoeff()
{
	size_t n_landmarks = pars_.landmark_src.size();
	tar_landmarks_.resize(n_landmarks*3); 
	RowMajorSparseMatrix landmark_coeff(n_landmarks*3, 12*num_sample_nodes);
	std::vector<Triplet> coeffs;
	for(int idx = 0; idx < n_landmarks; idx++)
	{
		int src_idx = pars_.landmark_src[idx];
		for (int k = 0; k < 3; k++)
		{
			for (RowMajorSparseMatrix::InnerIterator it(align_coeff_PV0_, src_idx*3+k); it; ++it)
			{
				coeffs.push_back(Triplet(idx*3+k, it.col(), it.value()));
			}
			tar_landmarks_[idx*3+k] = tar_points_(k, pars_.landmark_tar[idx]) - nodes_P_[src_idx*3+k];
		}
	}
	landmark_coeff.setFromTriplets(coeffs.begin(), coeffs.end());
	landmark_mul_ = landmark_coeff.transpose() * landmark_coeff; 
	landmark_right_ = landmark_coeff.transpose() * tar_landmarks_; 

	std::vector<Triplet> coeffs_fine;
	landmark_right_fine_.resize(n_src_vertex_*3);
	landmark_mul_fine_.resize(n_src_vertex_*3, n_src_vertex_*3);
	landmark_right_fine_.setZero();
	for(int idx = 0; idx < n_landmarks; idx++)
	{
		int src_idx = pars_.landmark_src[idx];
		for (int k = 0; k < 3; k++)
		{
			coeffs_fine.push_back(Triplet(src_idx*3+k, src_idx*3+k, 1.0));
			landmark_right_fine_[src_idx*3+k] = tar_points_(k, pars_.landmark_tar[idx]);
		}
	}
	landmark_mul_fine_.setFromTriplets(coeffs_fine.begin(), coeffs_fine.end());
}

void NonrigidSpareRegistration::InitWelschParam()
{
    // welsch parameters
    weight_d_.resize(n_src_vertex_*3);
    weight_d_.setOnes();

    // Initialize correspondences
    InitCorrespondence(correspondence_pairs_);

    VectorX init_nus(correspondence_pairs_.size());
#pragma omp parallel for
    for(size_t i = 0; i < correspondence_pairs_.size(); i++)
    {
        Vector3 closet = correspondence_pairs_[i].position;
        init_nus[i] = (src_mesh_->point(src_mesh_->vertex_handle(correspondence_pairs_[i].src_idx))
                       - Vec3(closet[0], closet[1], closet[2])).norm();
    }
    igl::median(init_nus, pars_.Data_nu);

    if(pars_.calc_gt_err&&n_src_vertex_ == n_tar_vertex_)
    {
        VectorX gt_err(n_src_vertex_);
        for(int i = 0; i < n_src_vertex_; i++)
        {
            gt_err[i] = (src_mesh_->point(src_mesh_->vertex_handle(i)) - tar_mesh_->point(tar_mesh_->vertex_handle(i))).norm();
        }
        pars_.init_gt_mean_errs = std::sqrt(gt_err.squaredNorm()/n_src_vertex_);
        pars_.init_gt_max_errs = gt_err.maxCoeff();
    }
}
void NonrigidSpareRegistration::InitwithLandmark()
{
	// Scalar energy=0., landmark_err=0., reg_err=0., rot_err=0., arap_err=0.;
	VectorX prevV = VectorX::Zero(n_src_vertex_ * 3);

	// welsch_sweight
 	bool run_once = true;

	Timer time;
	Timer::EventID begin_time, run_time;


	w_smo = optimize_w_smo;

	// pars_.each_energys.push_back(0.0);
	// pars_.each_gt_max_errs.push_back(pars_.init_gt_max_errs);
	// pars_.each_gt_mean_errs.push_back(pars_.init_gt_mean_errs);
	// pars_.each_iters.push_back(0);
	// pars_.each_times.push_back(pars_.non_rigid_init_time);
	// pars_.each_term_energy.push_back(Vector4(0, 0, 0, 0));

	Scalar gt_err;

	VectorX prev_X = X_;

	begin_time = time.get_time();
	
	#ifdef DEBUG
	double find_cp_time = 0.0;
	double construct_mat_time = 0.0;
	double solve_eq_time = 0.0;
	double calc_energy_time = 0.0;
	double robust_weight_time = 0.0;
	double update_r_time = 0.0;
	#endif

	std::cout << "init A" << std::endl;
	RowMajorSparseMatrix A_fixed_coeff = optimize_w_smo * reg_coeff_B_.transpose() * reg_cwise_weights_.asDiagonal() * reg_coeff_B_  + optimize_w_rot * rigid_coeff_L_ + optimize_w_landmark * landmark_mul_; // + optimize_w_arap * arap_coeff_mul_;

	int out_iter = 0;
	while (out_iter < pars_.max_outer_iters)
	{
		#ifdef DEBUG
		Timer::EventID inner_start_time = time.get_time();
		#endif

		// int welsch_iter;
		int total_inner_iters = 0;

		// std::cout << "out_iter" << out_iter << std::endl;

		mat_A0_ = A_fixed_coeff; 
		vec_b_ = optimize_w_smo * reg_coeff_B_.transpose() * reg_cwise_weights_.asDiagonal() * reg_right_D_ + optimize_w_rot * rigid_coeff_J_ * nodes_R_ + optimize_w_landmark * landmark_right_; // + optimize_w_arap * arap_coeff_.transpose() * arap_right_;

		#ifdef DEBUG
		Timer::EventID end_construct_eq = time.get_time();
		eps_time1 = time.elapsed_time(end_robust_weight, end_construct_eq);
		construct_mat_time += eps_time1;	
		#endif

		if (run_once)
		{
			solver_.analyzePattern(mat_A0_);
			run_once = false;
		}
		solver_.factorize(mat_A0_);
		X_ = solver_.solve(vec_b_);		

		run_time = time.get_time();
		double eps_time = time.elapsed_time(begin_time, run_time);
	
		#ifdef DEBUG
		eps_time1 = time.elapsed_time(end_construct_eq, run_time);
		solve_eq_time += eps_time1;
		#endif


		#ifdef DEBUG
		energy = CalcEnergy(align_err, reg_err, rot_err, arap_err, reg_cwise_weights_);
		// std::cout << "energy = " << energy << std::endl;
		#endif 
		// std::cout << "X = " << X_ << X_.sum() << std::endl;
		deformed_points_ = align_coeff_PV0_ * X_ + nodes_P_;
		// std::cout << "diff = " << (deformed_points_ - prevV).norm() << std::endl;
		
		#ifdef DEBUG
		Timer::EventID end_calc_energy = time.get_time();
		eps_time1 = time.elapsed_time(run_time, end_calc_energy);
		calc_energy_time += eps_time1;
		#endif

		

		#ifdef DEBUG
		Timer::EventID end_update_r = time.get_time();
		eps_time1 = time.elapsed_time(end_calc_energy, end_update_r);
		update_r_time += eps_time1;
		#endif
		
		CalcLocalRotations(0);
		CalcNodeRotations();
		CalcDeformedNormals();

		if (n_src_vertex_ == n_tar_vertex_)
			gt_err = (deformed_points_ - Eigen::Map<VectorX>(tar_points_.data(), 3 * n_src_vertex_)).squaredNorm();

		// save results
		// pars_.each_gt_mean_errs.push_back(gt_err);
		// pars_.each_gt_max_errs.push_back(0);
		// pars_.each_energys.push_back(energy);
		// pars_.each_iters.push_back(total_inner_iters);
		// pars_.each_term_energy.push_back(Vector4(0.0, reg_err, rot_err, arap_err));

		if((deformed_points_ - prevV).norm()/sqrtf(n_src_vertex_) < pars_.stop_coarse)
		{
			break;
		}
		prevV = deformed_points_;
		out_iter++;

		#ifdef DEBUG
		Timer::EventID end_find_cp2 = time.get_time();
		eps_time1 = time.elapsed_time(end_calc_energy, end_find_cp2);
		find_cp_time += eps_time1;
		#endif
	}

	#ifdef DEBUG
	std::cout << "find cp time = " << find_cp_time
	 << "\nconstruct_mat_timem = " << construct_mat_time
	 << "\nsolve_eq_time = " << solve_eq_time
	 << "\ncalc_energy_time = " << calc_energy_time
	 << "\nrobust_weight_time = " << robust_weight_time
	 << "\nupdate_r_time = " << update_r_time
	 << "\nacculate_iter = " << out_iter << std::endl;
	#endif
}


Scalar NonrigidSpareRegistration::DoNonRigid()
{
    // Data term parameters
    Scalar nu1 = pars_.Data_initk * pars_.Data_nu;
    if(pars_.landmark_src.size()>0)
	{
		optimize_w_landmark = 1.0/pars_.landmark_src.size();
		optimize_w_smo = pars_.w_smo/reg_coeff_B_.rows(); 
		optimize_w_rot = pars_.w_rot/num_sample_nodes;
		optimize_w_arap = pars_.w_arap_coarse/arap_coeff_.rows();
		std::cout << "Init Stage: optimize w_landmark = " << optimize_w_landmark << " | w_smo = " << optimize_w_smo << " | w_rot = " << optimize_w_rot  << " | w_arap = " << optimize_w_arap << std::endl;
		InitwithLandmark();
	}
    optimize_w_align = 1.0;
    //A parameter used to balance the influence of each energy component on the total energy.
    optimize_w_smo = pars_.w_smo/reg_coeff_B_.rows() * sampling_indices_.size();
    optimize_w_rot = pars_.w_rot/num_sample_nodes * sampling_indices_.size();
    optimize_w_arap = pars_.w_arap_coarse/arap_coeff_.rows() * sampling_indices_.size();
    std::cout << "Coarse Stage: optimize w_align = " << optimize_w_align << " | w_smo = " << optimize_w_smo << " | w_rot = " << optimize_w_rot << " | w_arap = " << optimize_w_arap << std::endl;
    GraphCoarseReg(nu1);//nu1 is the parameter of welsch kernel


    optimize_w_align = 1.0;
    optimize_w_arap = pars_.w_arap_fine /arap_coeff_fine_.rows() * n_src_vertex_;
    std::cout << "Fine Stage: optimize w_align = " << optimize_w_align << " | w_arap = " << optimize_w_arap << std::endl;
    PointwiseFineReg(nu1);//nu1 is still the parameter of welsch kernel


    // set the final deformed points to the source mesh
    //and return the final gt error if needed,if not return -1
    deformed_points_3X_ = Eigen::Map<Eigen::Matrix<double, 3, Eigen::Dynamic>>(deformed_points_.data(), 3, n_src_vertex_);
    return 0;
}


void NonrigidSpareRegistration::CalcNodeRotations()
{
#pragma omp parallel for
    for (int i = 0; i < num_sample_nodes; i++)
    {
        Matrix33 rot;
        Eigen::JacobiSVD<Matrix33> svd(Eigen::Map<Matrix33>(X_.data()+12*i, 3,3), Eigen::ComputeFullU | Eigen::ComputeFullV);
        if (svd.matrixU().determinant()*svd.matrixV().determinant() < 0.0) {
            Vector3 S = Vector3::Ones(); S(2) = -1.0;
            rot = svd.matrixU()*S.asDiagonal()*svd.matrixV().transpose();
        }
        else {
            rot = svd.matrixU()*svd.matrixV().transpose();
        }
        nodes_R_.segment(9 * i, 9) = Eigen::Map<VectorX>(rot.data(), 9);
    }
}




void NonrigidSpareRegistration::FullInARAPCoeff()
{
    arap_laplace_weights_.resize(n_src_vertex_);
    Timer timer;

    if(pars_.use_geodesic_dist)
    {
        for (int i = 0; i < n_src_vertex_; i++)
        {
            int nn = 0;
            OpenMesh::VertexHandle vh = src_mesh_->vertex_handle(i);
            for (auto vv = src_mesh_->vv_begin(vh); vv != src_mesh_->vv_end(vh); vv++)
            {
                nn++;
            }
            arap_laplace_weights_[i] = 1.0 / nn;
        }

        std::vector<Triplet> coeffs;
        std::vector<Triplet> coeffs_fine;
        for (int i = 0; i < src_mesh_->n_halfedges(); i++)
        {
            int src_idx = src_mesh_->from_vertex_handle(src_mesh_->halfedge_handle(i)).idx();
            int tar_idx = src_mesh_->to_vertex_handle(src_mesh_->halfedge_handle(i)).idx();
            Scalar w = sqrtf(arap_laplace_weights_[src_idx]);


            for (int k = 0; k < 3; k++)
            {
                for (RowMajorSparseMatrix::InnerIterator it(align_coeff_PV0_, src_idx*3+k); it; ++it)
                {
                    coeffs.push_back(Triplet(i*3+k, it.col(), w*it.value()));
                }
                for (RowMajorSparseMatrix::InnerIterator it(align_coeff_PV0_, tar_idx*3+k); it; ++it)
                {
                    coeffs.push_back(Triplet(i*3+k, it.col(), -w*it.value()));
                }
            }


            coeffs_fine.push_back(Triplet(i * 3, src_idx * 3, w));
            coeffs_fine.push_back(Triplet(i * 3 + 1, src_idx * 3 + 1, w));
            coeffs_fine.push_back(Triplet(i * 3 + 2, src_idx * 3 + 2, w));
            coeffs_fine.push_back(Triplet(i * 3, tar_idx * 3, -w));
            coeffs_fine.push_back(Triplet(i * 3 + 1, tar_idx * 3 + 1, -w));
            coeffs_fine.push_back(Triplet(i * 3 + 2, tar_idx * 3 + 2, -w));
        }


        arap_coeff_.resize(src_mesh_->n_halfedges()*3, num_sample_nodes * 12);
        arap_coeff_.setFromTriplets(coeffs.begin(), coeffs.end());
        arap_coeff_mul_ = arap_coeff_.transpose() * arap_coeff_;


        arap_coeff_fine_.resize(src_mesh_->n_halfedges() * 3, n_src_vertex_ * 3);
        arap_coeff_fine_.setFromTriplets(coeffs_fine.begin(), coeffs_fine.end());
        arap_coeff_mul_fine_ = arap_coeff_fine_.transpose() * arap_coeff_fine_;
    }
    else
    {
        int nn = src_knn_indices_.rows();
        for (int i = 0; i < n_src_vertex_; i++)
        {
            arap_laplace_weights_[i] = 1.0 / nn;
        }

        std::vector<Triplet> coeffs;
        std::vector<Triplet> coeffs_fine;
        for(int src_idx = 0; src_idx < n_src_vertex_; src_idx++)
        {
            for(int j = 0; j < nn; j++)
            {
                int i = src_idx*nn+j;
                int tar_idx = src_knn_indices_(j, src_idx);
                Scalar w = sqrtf(arap_laplace_weights_[src_idx]);


                for (int k = 0; k < 3; k++)
                {
                    for (RowMajorSparseMatrix::InnerIterator it(align_coeff_PV0_, src_idx*3+k); it; ++it)
                    {
                        coeffs.push_back(Triplet(i*3+k, it.col(), w*it.value()));
                    }
                    for (RowMajorSparseMatrix::InnerIterator it(align_coeff_PV0_, tar_idx*3+k); it; ++it)
                    {
                        coeffs.push_back(Triplet(i*3+k, it.col(), -w*it.value()));
                    }
                }


                coeffs_fine.push_back(Triplet(i * 3, src_idx * 3, w));
                coeffs_fine.push_back(Triplet(i * 3 + 1, src_idx * 3 + 1, w));
                coeffs_fine.push_back(Triplet(i * 3 + 2, src_idx * 3 + 2, w));
                coeffs_fine.push_back(Triplet(i * 3, tar_idx * 3, -w));
                coeffs_fine.push_back(Triplet(i * 3 + 1, tar_idx * 3 + 1, -w));
                coeffs_fine.push_back(Triplet(i * 3 + 2, tar_idx * 3 + 2, -w));
            }
        }

        arap_coeff_.resize(n_src_vertex_*nn*3, num_sample_nodes * 12);
        arap_coeff_.setFromTriplets(coeffs.begin(), coeffs.end());
        arap_coeff_mul_ = arap_coeff_.transpose() * arap_coeff_;


        arap_coeff_fine_.resize(n_src_vertex_*nn * 3, n_src_vertex_ * 3);
        arap_coeff_fine_.setFromTriplets(coeffs_fine.begin(), coeffs_fine.end());
        arap_coeff_mul_fine_ = arap_coeff_fine_.transpose() * arap_coeff_fine_;
    }
}


void NonrigidSpareRegistration::CalcARAPRight()
{
    if(pars_.use_geodesic_dist)
    {
#pragma omp parallel for
        for (int i = 0; i < src_mesh_->n_halfedges(); i++)
        {
            int src_idx = src_mesh_->from_vertex_handle(src_mesh_->halfedge_handle(i)).idx();
            int tar_idx = src_mesh_->to_vertex_handle(src_mesh_->halfedge_handle(i)).idx();

            Vector3 vij = local_rotations_.block(0, 3 * src_idx,3, 3) * (src_points_.col(src_idx) - src_points_.col(tar_idx));

            Scalar w = sqrtf(arap_laplace_weights_[src_idx]);


            arap_right_[i * 3] = w*(vij[0] - nodes_P_[src_idx * 3] + nodes_P_[tar_idx * 3]);
            arap_right_[i * 3 + 1] = w*(vij[1] - nodes_P_[src_idx * 3 + 1] + nodes_P_[tar_idx * 3 + 1]);
            arap_right_[i * 3 + 2] = w*(vij[2] - nodes_P_[src_idx * 3 + 2] + nodes_P_[tar_idx * 3 + 2]);
        }
    }
    else
    {
        int nn = src_knn_indices_.rows();
#pragma omp parallel for
        for (int src_idx = 0; src_idx < n_src_vertex_; src_idx++)
        {
            for (int j = 0; j < nn; j++)
            {
                int i = src_idx*nn + j;
                int tar_idx = src_knn_indices_(j, src_idx);

                Vector3 vij = local_rotations_.block(0, 3 * src_idx,3, 3) * (src_points_.col(src_idx) - src_points_.col(tar_idx));

                Scalar w = sqrtf(arap_laplace_weights_[src_idx]);


                arap_right_[i * 3] = w*(vij[0] - nodes_P_[src_idx * 3] + nodes_P_[tar_idx * 3]);
                arap_right_[i * 3 + 1] = w*(vij[1] - nodes_P_[src_idx * 3 + 1] + nodes_P_[tar_idx * 3 + 1]);
                arap_right_[i * 3 + 2] = w*(vij[2] - nodes_P_[src_idx * 3 + 2] + nodes_P_[tar_idx * 3 + 2]);
            }
        }
    }
}

void NonrigidSpareRegistration::CalcARAPRightFine()
{
    if(pars_.use_geodesic_dist)
    {
#pragma omp parallel for
        for (int i = 0; i < src_mesh_->n_halfedges(); i++)
        {
            int src_idx = src_mesh_->from_vertex_handle(src_mesh_->halfedge_handle(i)).idx();
            int tar_idx = src_mesh_->to_vertex_handle(src_mesh_->halfedge_handle(i)).idx();

            Vector3 vij = local_rotations_.block(0, 3 * src_idx, 3, 3) * (src_points_.col(src_idx) - src_points_.col(tar_idx));

            Scalar w = sqrtf(arap_laplace_weights_[src_idx]);

            arap_right_fine_[i * 3] = w*(vij[0]);
            arap_right_fine_[i * 3 + 1] = w*(vij[1]);
            arap_right_fine_[i * 3 + 2] = w*(vij[2]);
        }
    }
    else
    {
        int nn = src_knn_indices_.rows();
#pragma omp parallel for
        for(int src_idx = 0; src_idx < n_src_vertex_; src_idx++)
        {
            for(int j = 0; j < nn; j++)
            {
                int tar_idx = src_knn_indices_(j, src_idx);
                int i = src_idx*nn+j;

                Vector3 vij = local_rotations_.block(0, 3 * src_idx, 3, 3) * (src_points_.col(src_idx) - src_points_.col(tar_idx));

                Scalar w = sqrtf(arap_laplace_weights_[src_idx]);

                arap_right_fine_[i * 3] = w*(vij[0]);
                arap_right_fine_[i * 3 + 1] = w*(vij[1]);
                arap_right_fine_[i * 3 + 2] = w*(vij[2]);
            }
        }
    }
}


void NonrigidSpareRegistration::InitRotations()
{
    local_rotations_.resize(3, n_src_vertex_ * 3);
    local_rotations_.setZero();
#pragma omp parallel for
    for (int i = 0; i < n_src_vertex_; i++)
    {
        local_rotations_(0, i * 3) = 1;
        local_rotations_(1, i * 3 + 1) = 1;
        local_rotations_(2, i * 3 + 2) = 1;
    }
}


void NonrigidSpareRegistration::CalcLocalRotations(bool isCoarseAlign)
{

#pragma omp parallel for
    for (int i = 0; i < n_src_vertex_; i++)
    {
        Matrix33 sum;
        sum.setZero();

        int nn = 0;
        if(pars_.use_geodesic_dist)
        {
            OpenMesh::VertexHandle vh = src_mesh_->vertex_handle(i);
            for (auto vv = src_mesh_->vv_begin(vh); vv != src_mesh_->vv_end(vh); vv++)
            {
                int neighbor_idx = vv->idx();
                Vector3 dv = src_points_.col(i) - src_points_.col(neighbor_idx);
                Vector3 new_dv = deformed_points_.segment(3 * i, 3) - deformed_points_.segment(3 * neighbor_idx, 3);
                sum += dv * new_dv.transpose();
                nn++;
            }
        }
        else
        {
            nn = src_knn_indices_.rows();
            for(int j = 0; j < nn; j++)
            {
                int neighbor_idx = src_knn_indices_(j, i);
                Vector3 dv = src_points_.col(i) - src_points_.col(neighbor_idx);
                Vector3 new_dv = deformed_points_.segment(3 * i, 3) - deformed_points_.segment(3 * neighbor_idx, 3);
                sum += dv * new_dv.transpose();
            }
        }


        sum*= 1.0*optimize_w_arap/nn;

        if(!isCoarseAlign)
        {
            int tar_idx = correspondence_pairs_[i].tar_idx;
            Vector3 d = deformed_points_.segment(3 * i, 3) - tar_points_.col(tar_idx);
            Scalar c = (target_normals_.col(tar_idx) + deformed_normals_.col(i)).dot(d);
            Scalar d_norm2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
            Vector3 h = deformed_normals_.col(i) - c*d/d_norm2;

            Scalar w = optimize_w_align*d_norm2*weight_d_[i];
            sum += w * src_normals_.col(i) * h.transpose();
        }
        else if(vertex_sample_indices_[i] >= 0)
        {
            int tar_idx = correspondence_pairs_[vertex_sample_indices_[i]].tar_idx;
            Vector3 d = deformed_points_.segment(3 * i, 3) - tar_points_.col(tar_idx);
            Scalar c = (target_normals_.col(tar_idx) + deformed_normals_.col(i)).dot(d);
            Scalar d_norm2 = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
            Vector3 h = deformed_normals_.col(i) - c*d/d_norm2;

            Scalar w = optimize_w_align*d_norm2*weight_d_[vertex_sample_indices_[i]];
            sum += w * src_normals_.col(i) * h.transpose();
        }


        Eigen::JacobiSVD<Matrix33> svd(sum, Eigen::ComputeFullU | Eigen::ComputeFullV);

        if (svd.matrixU().determinant()*svd.matrixV().determinant() < 0.0) {
            Vector3 S = Vector3::Ones(); S(2) = -1.0;
            sum = svd.matrixV()*S.asDiagonal()*svd.matrixU().transpose();
        }
        else {
            sum = svd.matrixV()*svd.matrixU().transpose();
        }
        for (int s = 0; s < 3; s++)
        {
            for (int t = 0; t < 3; t++)
            {
                local_rotations_(s, 3 * i + t) = sum(s, t);
            }
        }
    }
}

void NonrigidSpareRegistration::CalcDeformedNormals()
{
#pragma omp parallel for
    for (int i = 0; i < n_src_vertex_; i++)
    {
        deformed_normals_.col(i) = local_rotations_.block(0, i * 3, 3, 3) * src_normals_.col(i);
    }
}

void NonrigidSpareRegistration::CalcNormalsSum()
{
#pragma omp parallel for
    for(int i = 0; i < correspondence_pairs_.size(); i++)
    {
        int sidx = correspondence_pairs_[i].src_idx;
        int tidx = correspondence_pairs_[i].tar_idx;
        int j = 0;
        for (RowMajorSparseMatrix::InnerIterator it(normals_sum_, i); it; ++it)
        {
            it.valueRef() = deformed_normals_(j, sidx) + target_normals_(j, tidx);
            j++;
        }
    }
}


void NonrigidSpareRegistration::InitNormalsSum()
{
    std::vector<Triplet> coeffs(3 * correspondence_pairs_.size());
    normals_sum_.resize(correspondence_pairs_.size(), 3*n_src_vertex_);
    normals_sum_.setZero();

#pragma omp parallel for
    for(int i = 0; i < correspondence_pairs_.size(); i++)
    {
        int sidx = correspondence_pairs_[i].src_idx;
        int tidx = correspondence_pairs_[i].tar_idx;
        coeffs[i * 3] = Triplet(i, 3 * sidx, deformed_normals_(0, sidx) + target_normals_(0, tidx));
        coeffs[i * 3 + 1] = Triplet(i, 3 * sidx + 1, deformed_normals_(1, sidx) + target_normals_(1, tidx));
        coeffs[i * 3 + 2] = Triplet(i, 3 * sidx + 2, deformed_normals_(2, sidx) + target_normals_(2, tidx));
    }
    normals_sum_.setFromTriplets(coeffs.begin(), coeffs.end());
}

void NonrigidSpareRegistration::PointwiseFineReg(Scalar nu1)
{

    Scalar energy=-1., align_err=-1., arap_err=-1.;

    VectorX prevV = VectorX::Zero(n_src_vertex_ * 3);

    bool run_once = true;


    // Smooth term parameters
    w_align = optimize_w_align;
    w_smo = optimize_w_smo;
    w_align = optimize_w_align *(2.0*nu1*nu1);


    Scalar gt_err = -1;

    double construct_mat_time = 0.0;
    double solve_eq_time = 0.0;


    int out_iter = 0;
    while (out_iter < pars_.max_outer_iters)
    {
        // Find clost points
        FindClosestPoints(target_tree_,deformed_points_,correspondence_pairs_);
        // according correspondence_pairs to update corres_U0_;
        corres_U0_.setZero();
        weight_d_.resize(n_src_vertex_);


        for (size_t i = 0; i < correspondence_pairs_.size(); i++)
        {
            corres_U0_.segment(i * 3, 3) = correspondence_pairs_[i].position;
            weight_d_[i] = correspondence_pairs_[i].min_dist2;
            int tar_idx = correspondence_pairs_[i].tar_idx;
            if(deformed_normals_.col(i).dot(target_normals_.col(tar_idx))<0)
                weight_d_[i] = -1;
        }

        // update weight

        welsch_weight(weight_d_, nu1);

        if (run_once == true && pars_.use_landmark == true)
        {
            weight_d_.setOnes();
        }

        // construct matrix A0 and pre-decompose
        if(out_iter==0)
            InitNormalsSum();
        else
            CalcNormalsSum();

        RowMajorSparseMatrix normals_sum_mul = normals_sum_.transpose() * weight_d_.asDiagonal()* normals_sum_;
        mat_A0_ = optimize_w_align * normals_sum_mul
                  + optimize_w_arap * arap_coeff_mul_fine_;

        CalcARAPRightFine();

        vec_b_ = optimize_w_align * normals_sum_mul * corres_U0_
                 + optimize_w_arap * arap_coeff_fine_.transpose() * arap_right_fine_;


        if (run_once)
        {
            solver_.analyzePattern(mat_A0_);
            run_once = false;
        }
        solver_.factorize(mat_A0_);
        deformed_points_ = solver_.solve(vec_b_);
        CalcLocalRotations(false);
        CalcDeformedNormals();
        if((deformed_points_ - prevV).norm()/sqrtf(n_src_vertex_) < pars_.stop_fine)
        {
            break;
        }
        prevV = deformed_points_;
        out_iter++;
    }
}

void NonrigidSpareRegistration::GraphCoarseReg(Scalar nu1)
{
    Scalar energy=0., align_err=0., reg_err=0., rot_err=0., arap_err=0.;
    VectorX prevV = VectorX::Zero(n_src_vertex_ * 3);

    // welsch_sweight
    bool run_once = true;


    w_align = optimize_w_align;
    w_smo = optimize_w_smo; 
    w_align = optimize_w_align *(2.0*nu1*nu1);

    Scalar gt_err;

    VectorX prev_X = X_;

    RowMajorSparseMatrix A_fixed_coeff = optimize_w_smo * reg_coeff_B_.transpose() * reg_cwise_weights_.asDiagonal() * reg_coeff_B_  + optimize_w_rot * rigid_coeff_L_ + optimize_w_arap * arap_coeff_mul_;
    int out_iter = 0;
    while (out_iter < pars_.max_outer_iters)
    {


        correspondence_pairs_.clear();
        FindClosestPoints(target_tree_,correspondence_pairs_, deformed_points_, sampling_indices_);
        corres_U0_.setZero();
        weight_d_.resize(correspondence_pairs_.size());
        weight_d_.setConstant(-1);

#pragma omp parallel for
        for (size_t i = 0; i < correspondence_pairs_.size(); i++)
        {
            corres_U0_.segment(correspondence_pairs_[i].src_idx * 3, 3) = correspondence_pairs_[i].position;
            weight_d_[i] = correspondence_pairs_[i].min_dist2;
            if(deformed_normals_.col(correspondence_pairs_[i].src_idx).dot(target_normals_.col(correspondence_pairs_[i].tar_idx))<0)
                weight_d_[i] = -1;
        }

        // update weight

        welsch_weight(weight_d_, nu1);

        if (run_once == true && pars_.use_landmark == true)
        {
            weight_d_.setOnes();
        }

        if(out_iter==0)
            InitNormalsSum();
        else
            CalcNormalsSum();

        diff_UP_ = (corres_U0_ - nodes_P_);

        RowMajorSparseMatrix weight_NPV = normals_sum_ * align_coeff_PV0_;

        mat_A0_ = optimize_w_align * weight_NPV.transpose() * weight_d_.asDiagonal() *  weight_NPV + A_fixed_coeff;

        CalcARAPRight();
        vec_b_ = optimize_w_align * weight_NPV.transpose() * weight_d_.asDiagonal() * normals_sum_ * diff_UP_ + optimize_w_smo * reg_coeff_B_.transpose() * reg_cwise_weights_.asDiagonal() * reg_right_D_ + optimize_w_rot * rigid_coeff_J_ * nodes_R_ + optimize_w_arap * arap_coeff_.transpose() * arap_right_;


        if (run_once)
        {
            solver_.analyzePattern(mat_A0_);
            run_once = false;
        }
        solver_.factorize(mat_A0_);
        X_ = solver_.solve(vec_b_);
        deformed_points_ = align_coeff_PV0_ * X_ + nodes_P_;
        CalcLocalRotations(true);
        CalcNodeRotations();
        CalcDeformedNormals();

        if (n_src_vertex_ == n_tar_vertex_)
            gt_err = (deformed_points_ - Eigen::Map<VectorX>(tar_points_.data(), 3 * n_src_vertex_)).squaredNorm();


        if((deformed_points_ - prevV).norm()/sqrtf(n_src_vertex_) < pars_.stop_coarse)
        {
            break;
        }
        prevV = deformed_points_;
        out_iter++;
    }

}

void NonrigidSpareRegistration::Paras_init(
    int iters,
    double stopcoarse,
    double stopfine,
    bool uselandmark,
    std::vector<int> src,
    std::vector<int> tar)
{
    paras.max_outer_iters = iters;
    paras.stop_coarse = stopcoarse;
    paras.stop_fine = stopfine;

    if (uselandmark && !src.empty() && src.size() == tar.size()) {
        paras.use_landmark = true;
        paras.landmark_src = src;
        paras.landmark_tar = tar;
    } else {
        paras.use_landmark = false;
        paras.landmark_src.clear();
        paras.landmark_tar.clear();
    }
}


void NonrigidSpareRegistration::Register()
{
	if(src_mesh.n_vertices()==0 || tar_mesh.n_vertices()==0)
        exit(0);

    if(src_mesh.n_vertices() != tar_mesh.n_vertices())
        paras.calc_gt_err = false;

    if(src_mesh.n_faces()==0)
        paras.use_geodesic_dist = false;

    if(paras.use_landmark)
        read_landmark(landmark_file.c_str(), paras.landmark_src, paras.landmark_tar);
    if(normalize)// The default value of `normalize` is `true`.
        scale = mesh_scaling(src_mesh, tar_mesh);
    pars_ = paras;
    Timer time;
    std::cout << "registration to initial... (mesh scale: " << scale << ")" << std::endl;
    Timer::EventID time1 = time.get_time();
    Init_data();
    Initialize();
    Timer::EventID time2 = time.get_time();
    std::cout << "non-rigid registration... (graph node number: " << pars_.num_sample_nodes << ")" << std::endl;
    DoNonRigid();
    Timer::EventID time3 = time.get_time();
    std::cout << "Registration done!\ninitialize time : "
              << time.elapsed_time(time1, time2) << " s \tnon-rigid reg running time = " << time.elapsed_time(time2, time3) << " s" << std::endl;
    return;
}

void NonrigidSpareRegistration::Reg(const std::string& file_target,
                       const std::string& file_source,
                       const std::string& out_path )
{
    Read_data(file_target, file_source);
    Register();
    Output_data(out_path, "spare");
}
