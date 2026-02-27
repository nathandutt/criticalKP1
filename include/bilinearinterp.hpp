#pragma once
#include <stdio.h>
#include <iostream>
#include <array>
#include <cmath>
#include <vector>
#include "point.hpp"

//Header to find zeros in a square
//with bilinear interpolation
using array = std::array<double, 4>;
using sol_array = std::array<double, 2>;

constexpr double det(const double a, const double b, const double c, const double d) noexcept{
    return a*b-c*d;
}

inline int roots(const double c0, const double c1, const double c2, sol_array& y_sols){
     double disc = c1*c1-4.*c0*c2;
     if(disc < 0.) return 0;
     double sol_1 = 1./(2.*c2) * (-c1 + std::sqrt(disc));
     double sol_2 = 1./(2.*c2) * (-c1 - std::sqrt(disc));
     int returned = 0;
     if((0.<=sol_1) && (sol_1 <= 1.)){
         y_sols[0] = sol_1; returned++;
     }
     if((0.<=sol_2) && (sol_2 <= 1.)){
         y_sols[1] = sol_2; returned++;
     }
//     std::cout << "got 2 roots!" << std::endl;
//     std::cout << sol_1 << " " << sol_2 << std::endl;
     return returned;
}

inline int addif(int idx, const double y_sol, const array& ga, sol_array& x_sols){
    //Adds x_root corresponding to y root, only if within correct range
    //returns 0 if none added, 1 if added
    double x_candidate = -1.*(ga[0]+y_sol*ga[2])/(ga[1] + y_sol*ga[3]);
    if((x_candidate < 0.) || (x_candidate > 1.))
        return 0;
    x_sols[idx] = x_candidate;
    return 1;
}

//Supposes you've rescaled, so corners are (0, 0), (1, 0), (1, 1), (0,1)
inline int getReducedZeros(const array& ga, const array& gb, sol_array& x_sols, sol_array& y_sols){
    //Finds zero of bilinear interpolation.
    //Stores them in sols, returns number of solutions;
    double c0 = det(ga[0], gb[1], ga[1], gb[0]);
    double c1 = det(ga[0], gb[3], ga[3], gb[0]) + det(ga[2], gb[1], ga[1], gb[2]);
    double c2 = det(ga[2], gb[3], ga[3], gb[2]);

    int s = roots(c0, c1, c2, y_sols);
    int returned = 0;
    switch(s){
        case 0:
            return 0;
        case 1:
            returned += addif(0, y_sols[0], ga, x_sols);
            break;
        case 2:
            returned += addif(0, y_sols[0], ga, x_sols);
            returned += addif(1, y_sols[1], ga, x_sols);
            break;
    }
    return returned;
}

//Now, imagine we have 4 corners, and the 4 values
//Supposes array correspond to points (0, 0), (0, 1), (1, 0), (1,1) in terms of corners
inline std::vector<Point> getZeros(const sol_array& xc, const sol_array& yc, const array& phi_x, const array& phi_y){
    array gx{}; array gy{};
    //Construct bilinear coefficient arrays
    gx[0] = phi_x[0];
    gx[1] = phi_x[1] - phi_x[0];
    gx[2] = phi_x[2] - phi_x[0];
    gx[3] = phi_x[3] + phi_x[0] - phi_x[1] - phi_x[2];

    gy[0] = phi_y[0];
    gy[1] = phi_y[1] - phi_y[0];
    gy[2] = phi_y[2] - phi_y[0];
    gy[3] = phi_y[3] + phi_y[0] - phi_y[1] - phi_y[2];

    sol_array x_s{}; sol_array y_s{};
    int num_zeros = getReducedZeros(gx, gy, x_s, y_s);

    //Zeros are in reduced form, we need to get them back
    for(int i = 0; i < num_zeros; i++){
        x_s[i] = x_s[i] * xc[0] + (1. - x_s[i])*xc[1];
        y_s[i] = y_s[i] * yc[0] + (1. - y_s[i])*yc[1];
    }

    std::vector<Point> res;
    for(int i=0; i < num_zeros; i++){
        res.emplace_back(Point(x_s[i], y_s[i]));
    }

    return res;
}
