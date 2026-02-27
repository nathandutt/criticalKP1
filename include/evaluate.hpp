#pragma once
#include "global.hpp"
#include <Eigen/Dense>

//
// Quickly computing x, y derivatives of field when given poles at certain y
// Computes over x, in a given range
Eigen::VectorXd Phi_x(const SavedState& state, const double x_i, const double x_f, const int points) {
    Eigen::VectorXd x_vals = Eigen::VectorXd::LinSpaced(points, x_i, x_f);
    Eigen::VectorXd f_vals = Eigen::VectorXd::Zero(points);

    for (int i = 0; i < state.poles.cols(); i++) {
        double a = state.poles(0, i);
        double b = state.poles(1, i);
        auto denom = ((x_vals.array() - a).square() + b * b).square(); // denominator
        f_vals += (((x_vals.array() - a).square() - b * b) / denom).matrix();
    }

    f_vals *= -1.;
    return f_vals;
}

Eigen::VectorXd Phi_y(const SavedState& state, const double x_i, const double x_f, const int points) {
    Eigen::VectorXd x_vals = Eigen::VectorXd::LinSpaced(points, x_i, x_f);
    Eigen::VectorXd f_vals = Eigen::VectorXd::Zero(points);

    for (int i = 0; i < state.poles.cols(); i++) {
        double a = state.poles(0, i);
        double b = state.poles(1, i);
        double k_re = state.velocity(0, i);
        double k_im = state.velocity(1, i);

        auto denom = ((x_vals.array() - a).square() + b * b).square();
        f_vals += (k_re * ((x_vals.array() - a).square() - b * b) / denom).matrix();
        f_vals += (k_im * 2.0 * (x_vals.array() - a) * b / denom).matrix();
    }

    return f_vals;
}
