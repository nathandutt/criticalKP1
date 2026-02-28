#include "polestate.hpp"
#include <stdexcept>
#include <iostream>

//Force computation
Column F(const Column& m, double& max_c) {
    // m.row(0) = x coordinates, m.row(1) = y coordinates
    Column f = Column::Zero();

    // Compute all pairwise differences
    Eigen::Array<double, N, N> dx = m.row(0).replicate(N,1).transpose() - m.row(0).replicate(N,1);
    Eigen::Array<double, N, N> dy = m.row(1).replicate(N,1).transpose() - m.row(1).replicate(N,1);

    // r_ij^2 = dx^2 + dy^2, avoid diagonal (i==j) by setting it to 1 temporarily
    Eigen::Array<double, N, N> r2 = dx.square() + dy.square();
    Eigen::Array<bool, N, N> diag_mask = (r2 == 0);
    r2 = r2 + diag_mask.cast<double>(); // temporarily prevent div by zero

    // Update max_c for adaptive timestep

    Eigen::Array<double, N, N> r6_inv = r2.pow(-3.0); // 1 / r^6

    // Numerators
    Eigen::Array<double, N, N> num_x = dx.cube() - 3*dx*dy.square();
    Eigen::Array<double, N, N> num_y = dy.cube() - 3*dy*dx.square();

    // Multiply by 1/r^6
    Eigen::Array<double, N, N> contrib_x = r6_inv * num_x;
    Eigen::Array<double, N, N> contrib_y = r6_inv * num_y;

    // Zero out diagonal (i==j)
    contrib_x = contrib_x * (!diag_mask).cast<double>();
    contrib_y = contrib_y * (!diag_mask).cast<double>();

    // Sum contributions for each i
    f.row(0) = contrib_x.colwise().sum();
    f.row(1) = contrib_y.colwise().sum();
    if(max_c < 0.){ 
	Eigen::Array<double, 1, N> magnitudes =
	    (f.row(0).array().square() + f.row(1).array().square()).sqrt();

	double max_force = magnitudes.maxCoeff();
	max_c = max_force;
    }
    f*=-2.;
    return f;
}

constexpr double adaptative_pow = 1.;
void PoleState::evolveRK4(const double timestep){
    //RK4 evolution with adaptative timestep
    double max_c = -1.;

    //first k
    Column k1_p = velocity;
    Column k1_v = F(poles, max_c);
    //Define adaptative timestep
    if(max_c > 1e6) 
	throw std::runtime_error("Numerical error");
    double astep = (max_c < 1.) ? timestep : timestep*std::pow(max_c, -adaptative_pow);

    Column k2_p = velocity + (0.5*astep)*k1_v;
    Column k2_v = F(poles + (0.5*astep)*k1_p, max_c);

    Column k3_p = velocity + (0.5*astep)*k2_v;
    Column k3_v = F(poles+(0.5*astep)*k2_p, max_c);

    Column k4_p = velocity + astep*k3_v;
    Column k4_v = F(poles + astep*k3_p, max_c);

    velocity += astep/6. * (k1_v + 2.*k2_v + 2.*k3_v + k4_v);
    poles += astep/6. * (k1_p + 2.*k2_p + 2.*k3_p + k4_p);
    y += astep;
}
