
#ifndef FRICP_H
#define FRICP_H

#include <AndersonAcceleration.h>
#define TUPLE_SCALE	  0.95
#define TUPLE_MAX_CNT 1000
#include "registration.h"






class RigidFricpRegistration:public Registration
{
public:
    RigidFricpRegistration(){};
    ~RigidFricpRegistration(){};
    std::string file_init_  ;
    VectorN source_mean_, target_mean_;
    double time;
    ICP::Parameters pars;
    void use_init_transform(bool a);
    AffineMatrix3 LogMatrix(const AffineMatrix3& T);
    Affine3d point_to_point(Matrix3X& X,
                         Matrix3X& Y,
                        const VectorX& w) ;
    Eigen::Affine3d point_to_plane(Eigen::Matrix3Xd& X,
                                   Eigen::Matrix3Xd& Y,
                                   const Eigen::Matrix3Xd& Norm,
                                   const Eigen::VectorXd& w,
                                   const Eigen::VectorXd& u) ;
    void point_to_point(Matrix3X& X, Matrix3X& Y,Matrix3X& Z, Vector3& source_mean_,
                        Vector3& target_mean_, ICP::Parameters& par);
    void point_to_plane(Eigen::Matrix3Xd& X,
                        Eigen::Matrix3Xd& Y, Eigen::Matrix3Xd& norm_x, Eigen::Matrix3Xd& norm_y,
                        Eigen::Vector3d& source_mean_, Eigen::Vector3d& target_mean_,
                        ICP::Parameters &par) ;
    void Paras_init(bool useinit =false,std::string fileinit=" ",int maxiter=100,double stop=1e-5);
    void Register()override;
    void Reg(const std::string& file_target,
                       const std::string& file_source,
                       const std::string& out_path)override;
    };


#endif
