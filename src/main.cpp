#include "global.hpp"
#include <algorithm>
#include "polestate.hpp"
#include "initialcondition.hpp"
#include "bilinearinterp.hpp"
#include "evaluate.hpp"
#include "npy.hpp"

#include <fstream>
#include <iostream>
#include <vector>
#include <complex>
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <iomanip>

using namespace std;
/*
 *       Format of params.txt
 *       t_i  t_f  t_step
 *       y_i  y_f  y_step y_write_step
 *       x_i  x_f  x_step 
 *       k1_re  k1_im  o1_re  o1_im
 *       ...
 *       kn_re  kn_im  on_re  on_im
 */

//Read Evolution parameters
struct Config{
    double y_i; double y_f; double y_step; double y_write_step;
    double x_i; double x_f; double x_points;
    double t_i; double t_f; double t_step;

    Config(fstream& file){
	file >> t_i >> t_f >> t_step;
	file >> y_i >> y_f >> y_step >> y_write_step;
	file >> x_i >> x_f >> x_points;
    }
    void Print(){
	cout << "Y parameters are: " << endl
	    << "Initial y: " << y_i << " Final y: " << y_f <<
	    " Pole integration step: " << y_step << " Zero finding step: " << y_write_step << endl
	    << "X parameters are: " << endl
	    << "Initial x: " << x_i << " Final x: " << x_f << " Zero finding step: " <<x_points << endl
	    << "T parameters are: " << endl
	    << "Initial t: " << t_i << " Final t: " << t_f << " Timestep:" << t_step << endl;
    }
};

//Read soliton parameters
void ReadParameters(fstream& file, vector<complex<double>>& k_s, vector<complex<double>>& offsets){
    double k_re, o_re, k_im, o_im;
    while(file >> k_re >> k_im >> o_re >> o_im){
	auto k = complex<double>(k_re, k_im); 
	auto o = complex<double>(o_re, o_im); 
	k_s.emplace_back(k);
	offsets.emplace_back(o);
	offsets.emplace_back(conj(o));
	k_s.emplace_back(conj(k));
    }
}

//Evolve soliton offsets for a fixed time
vector<complex<double>> EvolveOffsets(double t,const vector<complex<double>>& k_s, const vector<complex<double>>& offsets){
    auto new_offsets = std::vector<complex<double>>{};
    if(k_s.size()!=offsets.size()) throw std::runtime_error("Nonequal k_s and offsets");
    for(unsigned int i = 0; i < k_s.size(); i++){
	auto new_o = offsets[i] - 12. * k_s[i]*k_s[i]*t;
	new_offsets.emplace_back(new_o);
    }
    return new_offsets;
}

constexpr unsigned long max_points=10;
int addPointsToData(const std::vector<Point>& cpoints, std::vector<double>& data){
    //We add the critical points to data.
    //If there are less than 6, (To be adapted maybe later on for larger N), we add a padding
    //If more than six, just skip.
    if(cpoints.size() > max_points){
	cout << "Too many points for writing to numpy" << endl;
	return 0;
    }
    int emplaced = 0;

    for(const auto& p : cpoints){
	data.emplace_back(p.x); data.emplace_back(p.y);
	emplaced++;
    }
    //Now add dummy points
    for(int i = 0; i + emplaced < max_points; i++){
	data.emplace_back(std::numeric_limits<double>::quiet_NaN());
	data.emplace_back(std::numeric_limits<double>::quiet_NaN());
    }
    return 1;
}


// !!! We suppose params.txt is well formatted, see above
int main(){
    const string input_file = "params.txt";
    const string output_dir = "Output/";
    fstream file;
    auto k_s = vector<complex<double>>{};
    auto offsets = vector<complex<double>>{};

    file.open(input_file, ios::in);
    //Get config
    auto config = Config(file);

    cout << "Take a moment to check parameters please. " << endl;
    config.Print();

    //Get soliton parameters
    ReadParameters(file, k_s, offsets);
    for(const auto & k : k_s)
	cout << k << ", ";
    cout << endl;
    file.close();

    int pole_number = k_s.size();
    assert(pole_number == N);

    double current_time = config.t_i;
    double final_time = config.t_f;
    int steps = (final_time - current_time)/config.t_step + 1;
    
    std::vector<double> data;
    data.reserve(12*steps+5); // A priori two coordinates, and 6 critical points per step
    unsigned long written = 0;
    while(current_time < config.t_f){

	cout << "Time = " << current_time << endl;
	std::cout << "Evolving offsets" << endl;
	auto t_offsets = EvolveOffsets(current_time, k_s, offsets);
	cout << "Fetching initial condition" << endl;
	//Initialize poles
	auto [p, v] = InitialConditions(k_s, t_offsets, config.y_i);
	PoleState poles(config.y_i, p, v);
	for(int i = 0; i < N; i++)
	    cout << "("<<v(0,i) <<"," << v(1, i) <<")";
	cout << endl;
	//Vector to save
	int total_pt_estimate = floor((config.y_f - config.y_i)/config.y_write_step);
	std::vector<std::unique_ptr<SavedState>> states;
	states.reserve(1.5*total_pt_estimate);

	//Initialize loop parameters

	double current_y = config.y_i;
	double next_write_y = config.y_i;
	int write_step = config.y_write_step/config.y_step;
	int curr_step = 0;
	cout << "Starting y loop" << endl;
	Eigen::VectorXd x_vals = Eigen::VectorXd::LinSpaced(config.x_points, config.x_i, config.x_f);
	while(current_y < config.y_f){
	    if(curr_step%write_step ==0){
		poles.Insert(states);
		next_write_y += config.y_write_step;
	    }
	    curr_step++;
	    poles.Evolve(config.y_step);
	    current_y = poles.y;
	}
	cout << "Finished CMS evolution." << endl;

	//Finished evolution, now find critical points
	std::vector<Point> critical_points;
	critical_points.reserve(6);

	Eigen::VectorXd prev_phi_x, prev_phi_y;
	cout << "Now looking for critical pts" << endl; 
	if (!states.empty()) {
	    prev_phi_x = Phi_x(*states[0], config.x_i, config.x_f, config.x_points);
	    prev_phi_y = Phi_y(*states[0], config.x_i, config.x_f, config.x_points);
	}

	for (std::size_t i = 0; i + 1 < states.size(); ++i) {
	    SavedState& current = *states[i];
	    SavedState& next    = *states[i + 1];

	    // reuse prev_phi_x/prev_phi_y for current
	    Eigen::VectorXd phi_x_next = Phi_x(next, config.x_i, config.x_f, config.x_points);
	    Eigen::VectorXd phi_y_next = Phi_y(next, config.x_i, config.x_f, config.x_points);

	    double y1 = current.y;
	    double y2 = next.y;

	    for(int j = 0; j + 1 < config.x_points; j++){
		double maxx = std::max({phi_x_next(j), phi_x_next(j+1), prev_phi_x(j), prev_phi_x(j+1)});
		double minx = std::min({phi_x_next(j), phi_x_next(j+1), prev_phi_x(j), prev_phi_x(j+1)});
		double maxy = std::max({phi_y_next(j), phi_y_next(j+1), prev_phi_y(j), prev_phi_y(j+1)});
		double miny = std::min({phi_y_next(j), phi_y_next(j+1), prev_phi_y(j), prev_phi_y(j+1)});
		if( (maxx * minx) > 0. || ((miny * maxy) > 0.)) continue;
		std::array<double, 2> xc{x_vals(j), x_vals(j+1)};
		std::array<double, 2> yc{y1, y2};
		std::array<double, 4> phix{prev_phi_x(j), prev_phi_x(j+1), phi_x_next(j), phi_x_next(j+1)};
		std::array<double, 4> phiy{prev_phi_y(j), prev_phi_y(j+1), phi_y_next(j), phi_y_next(j+1)};
		std::vector<Point> new_points = getZeros(xc, yc, phix, phiy);
		for(const auto& pt : new_points)
		    critical_points.emplace_back(pt);
	    }
	    // Shift next to prev for next iteration
	    prev_phi_x = std::move(phi_x_next);
	    prev_phi_y = std::move(phi_y_next);
	}

	cout << "Finished. Got " << critical_points.size() << " critical points.";
	written += addPointsToData(critical_points, data);
	current_time += config.t_step;
    }

    //Writing to .npy
    string write_path = output_dir + "cpoints.npy";
    std::vector<unsigned long> shape{written, max_points, 2};
    const npy::npy_data_ptr<double> data_ptr{data.data(), shape, false};
    write_npy(write_path, data_ptr);
}
