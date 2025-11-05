#ifndef ROBUST_NORM_H_
#define ROBUST_NORM_H_

#include "registration.h" 

inline void welsch_weight(Eigen::VectorXd& r, double p)
{
    #pragma omp parallel for
    for (int i = 0; i < r.size(); ++i)
    {
        r[i] = (r[i] >= 0) ? std::exp(-r[i] / (2 * p * p)) : 0.0;
    }
}

inline double welsch_energy(const Eigen::VectorXd& r, double p)
{
    double energy = 0.0;
    #pragma omp parallel for reduction(+:energy)
    for (int i = 0; i < r.size(); ++i)
    {
        energy += 1.0 - std::exp(-r[i]  / (2 * p * p));
    }
    return energy;
}

inline void robust_weight(Eigen::VectorXd& r, double p)
{
    welsch_weight(r, p);
}

inline double get_energy(const Eigen::VectorXd& r, double p)
{
    return welsch_energy(r, p);
}


#endif
