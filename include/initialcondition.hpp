#pragma once
#include <complex>
#include <vector>
#include "polestate.hpp"
auto InitialConditions(const std::vector<std::complex<double>>& k_s, const std::vector<std::complex<double>>& offsets, const double y_i)
    ->std::pair<Column, Column>;


