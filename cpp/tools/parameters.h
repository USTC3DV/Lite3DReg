#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "Types.h"

//Parameter repository for different methods

// parameters for FRICP
namespace ICP{
    struct Parameters {
    double p;       /// paramter of the robust function/// para k
    int max_icp;    /// max ICP iteration
    double stop;    /// stopping criteria
    std::string out_path;
    int anderson_m;
    MatrixXX init_trans;
    MatrixXX gt_trans;
    bool has_groundtruth;
    double convergence_energy;
    double convergence_gt_mse;
    MatrixXX res_trans;
    double nu_begin_k;
    double nu_end_k;
    bool use_init;
    double nu_alpha;
    Parameters() :
        p(0.1),
        max_icp(100),
        stop(1e-5),
        anderson_m(5),
        has_groundtruth(false),
        convergence_energy(0.0),
        convergence_gt_mse(0.0),
        nu_begin_k(3),
        nu_end_k(1.0 / (3.0 * sqrt(3.0))),
        use_init(false),
        nu_alpha(1.0 / 2)
    {
        gt_trans = Eigen::Matrix4d::Identity();
        init_trans = Eigen::MatrixXd(); // 或指定大小
        res_trans = Eigen::MatrixXd();
    }
    };
}

namespace spare{
// parameters for spare non-rigid registration
struct Parameters
{
    int		max_outer_iters;    // nonrigid max iters
    Scalar	w_smo;              // smoothness weight 
    Scalar	w_rot;               // rotation matrix weight 
	Scalar  w_arap_coarse;      // ARAP weight for coarse alignment
    Scalar  w_arap_fine;        // ARAP weight for fine alignment 
    bool	use_landmark;
    bool    calc_gt_err;         // calculate ground truth error (DEBUG)
    std::vector<int> landmark_src;
    std::vector<int> landmark_tar;
    Scalar  Data_nu;
    Scalar  Data_initk;
    Scalar  stop_coarse;
    Scalar  stop_fine;
	// Sample para
    Scalar  uni_sample_radio;       // uniform sample radio
    bool    use_geodesic_dist;
    int     num_sample_nodes;
    // record the initial error
    Scalar  init_gt_mean_errs;
    Scalar  init_gt_max_errs;
    Scalar  mesh_scale;
    Parameters() // default
    {
        max_outer_iters = 30;
        w_smo = 0.01;  // smooth
        w_rot = 1e-4;   // orth
		w_arap_coarse = 500; // 10;
        w_arap_fine = 200;
        use_landmark = false;
        calc_gt_err = true;
        Data_nu = 0.0;
        Data_initk = 1;
        stop_coarse = 1e-3;
        stop_fine = 1e-4;
        // Sample para
        uni_sample_radio = 10;
        use_geodesic_dist = true;
        init_gt_mean_errs = .0;
    }
};
};

#endif // PARAMETERS_H