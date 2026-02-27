#pragma once
#include <memory>
#include "global.hpp"  

struct PoleState {
    double y;
    Column poles;
    Column velocity;

    
    PoleState() : poles(Column::Zero()), velocity(Column::Zero()) {}

    
    PoleState(const double y_i, const Column& p, const Column& v)
        : y(y_i), poles(p), velocity(v) {}

    
    PoleState(const double y_i, Column&& p, Column&& v)
        : y(y_i), poles(std::move(p)), velocity(std::move(v)) {}

    void Evolve(const double timestep);
    //implemented in .cpp, RK4 evolution of CalogeroMoser Hamiltonian

    void Insert(std::vector<std::unique_ptr<SavedState>>& v){
	v.emplace_back(std::make_unique<SavedState>(y, poles, velocity));
    }
     
};

