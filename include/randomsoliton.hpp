#pragma once
#include "global.hpp"
#include <vector>
#include <complex>
#include <random>

using vect = std::vector<std::complex<double>>;
using myRNG = std::mt19937;

//Header to generate random solitons, change N in global to 2*number of solitons wanted.

//A priori we must be careful that solitons don't get so close that CM integration fails or gives bad values
inline std::complex<double> singleRandomOffset(myRNG& rng, const double x_min, double x_max, double y_min, double y_max){
    std::uniform_real_distribution<> dis(0., 1.);
    double r1 = dis(rng);
    double r2 = dis(rng);
    double x = x_min*(1.-r1) + x_max*r1;
    double y = y_min*(1.-r2) + y_max*r2;
    return std::complex<double>(x, y);
}
inline vect randomOffsets(myRNG& rng, const double x_min, const double x_max, const double y_min, const double y_max){
    assert(N%2==0);
    vect offsets{};
    offsets.reserve(N);
    for(int i = 0; i<(N/2); i++){
	std::complex<double> offset = singleRandomOffset(rng, x_min, x_max, y_min, y_max);
	offsets.emplace_back(offset);
	offsets.emplace_back(std::conj(offset));
    }
    return offsets;
}

//Generate N//2 random k_s, only in upper half plane
//Generate with uniform distribution, if outside of circle just throw again
//Not sure if taking fabs hurts randomness, but not a very big concern.
inline std::complex<double> singleRandomK(myRNG& rng, const double max_mod){
    std::uniform_real_distribution<> dis(-max_mod, max_mod);
    double k_re = std::fabs(dis(rng));
    double k_im = dis(rng);
    if((k_re*k_re + k_im * k_im) > max_mod*max_mod) //We just try again
	return singleRandomK(rng, max_mod);
    return std::complex<double>(k_re, k_im);
}
inline vect randomKs(myRNG& rng, const double max_mod){
    assert(N%2==0);
    vect k_s{};
    k_s.reserve(N);
    for(int i = 0; i<(N/2); i++){
	std::complex<double> k = singleRandomK(rng, max_mod);
	k_s.emplace_back(k);
	k_s.emplace_back(std::conj(k));
    }
    return k_s;
}

