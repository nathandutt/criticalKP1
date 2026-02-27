#include <stdexcept>
#include "initialcondition.hpp"
#include "global.hpp"
using complex = std::complex<double>;
using Matrix  = Eigen::MatrixXcd;
using Vector  = Eigen::VectorXcd;

Matrix Hirota(const std::vector<complex>& k_s,
              const std::vector<complex>& offsets,
              const double y_i)
{
    if((int) k_s.size() != N || offsets.size() != N)
	throw std::runtime_error("Wrong number of k or offsets");
    Matrix H(N, N);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j)
                H(i, j) = 2.0 * k_s[i] * y_i + offsets[i];
            else
                H(i, j) = complex(0.0, 1.0) / (k_s[i] - k_s[j]);
        }
    }

    return H;
}

Matrix K(const std::vector<complex>& k_s)
{
    if((int) k_s.size() != N) 
	throw std::runtime_error("Wrong number of k or offsets");

    Matrix Km = Matrix::Zero(N, N);

    for (int i = 0; i < N; ++i)
        Km(i, i) = 2.0 * k_s[i];

    return Km;
}

auto InitialConditions(const std::vector<complex>& k_s,
                      const std::vector<complex>& offsets,
                      const double y_i)
-> std::pair<Column, Column>
{
    //Given soliton parameters at a fixed t. Compute inital pole positions and velocity, at a y value y_i.
    //This let's us not have to take asympotic values everytime, and we can start simulations at y=20, for example, where
    //Some deviation might already have happened

    if((int) k_s.size() != N || offsets.size() != N)
	throw std::runtime_error("Wrong number of k or offsets");

    Matrix H  = Hirota(k_s, offsets, y_i);
    Matrix Km = K(k_s);

    // Compute eigenvalues and right eigenvectors
    Eigen::ComplexEigenSolver<Matrix> solver(H);

    Vector eigenvalues_vec = solver.eigenvalues();
    Matrix VR              = solver.eigenvectors();

    // Compute zdot = diag(V^{-1} K V)
    Matrix Vinv = VR.inverse();
    Matrix res  = Vinv * Km * VR;

    Column pole;
    pole.setZero();
    Column velocity;
    velocity.setZero();

    for(int i = 0; i<N; ++i){
	pole(0, i) = eigenvalues_vec(i).real(); pole(1, i) = eigenvalues_vec(i).imag();
	velocity(0, i) = res(i, i).real(); velocity(1, i) = res(i, i).imag();
    }
    return {pole, velocity};
}
