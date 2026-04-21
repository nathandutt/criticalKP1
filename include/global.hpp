#pragma once
#include <Eigen/Dense>

constexpr int N = 4;
using Column = Eigen::Matrix<double, 2, N>;  

struct SavedState{
    double y;
    Column poles;
    Column velocity;

    SavedState(const double y, const Column& p, const Column& v) :
	y(y), poles(p), velocity(v) {}
};

