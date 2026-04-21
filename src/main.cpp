#include <sstream>
#include <fstream>
#include <iostream>

#include <iomanip>
#include <stdexcept>

#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>

#include "global.hpp" 			//Global parameters, set number of solitons here
					//Has to be a global constant <100 because uses stack memory

#include "randomsoliton.hpp"		//Generate random soliton parameters

#include "initialcondition.hpp"		//Get initial conditions for fixed t
					//(Diagonalization of Hirota matrix)

#include "polestate.hpp"		//RK4 evolution of poles

#include "bilinearinterp.hpp"		//Header defining interpolation functions for critical pts

#include "evaluate.hpp"			//Quick evaluation of field given poles

#include "npy.hpp"			//Header to easily write data to .npy format


using namespace std;
using myRNG = std::mt19937;

/*
 *       Format of params.txt
 *       t_i  t_end  t_step
 *       y_start  y_f  y_step y_write_step
 *       x_start  x_end  x_step 
 *       k1_re  k1_im  o1_re  o1_im
 *       ...
 *       kn_re  kn_im  on_re  on_im
 */

//Read Evolution parameters
struct Config{
    double y_start; double y_end; int y_points; double cms_step;
    double x_start; double x_end; int x_points;
    double t_start; double t_end; int t_points;

    Config(const std::string& filename){
	fstream file;
	file.open(filename, ios::in);

	string str;
	while(getline(file, str)){
	    stringstream ss(str);
	    
	    //Skip comments or empty lines
	    if(str[0] == '#' || str.empty())
		continue;
	    string key, val;
	    ss >> key >> val;

	    if(key == "y_start") y_start = stod(val);
	    else if(key == "y_end") y_end = stod(val);
	    else if(key == "y_points") y_points = stoi(val);
	    else if(key == "cms_step") cms_step = stod(val);
	    else if(key == "x_start") x_start = stod(val);
	    else if(key == "x_end") x_end = stod(val);
	    else if(key == "x_points") x_points = stoi(val);
	    else if(key == "t_start") t_start = stod(val);
	    else if(key == "t_end") t_end = stod(val);
	    else if(key == "t_points") t_points = stoi(val);
	    else{
		string exception = "Invalid Parameter name : " + key;
		throw std::runtime_error(exception);
	    }

	}
    }
    void Print(){
	cout << "Y parameters are: " << endl
	    << "Initial y: " << y_start << " Final y: " << y_end <<
	    " Pole integration step: " << cms_step << " Written y points: " << y_points << endl
	    << "X parameters are: " << endl
	    << "Initial x: " << x_start << " Final x: " << x_end << " x points :" <<x_points << endl
	    << "T parameters are: " << endl
	    << "Initial t: " << t_start << " Final t: " << t_end << "T points : " << t_points << endl;
    }
};

//Read soliton parameters
void readParameters(const std::string& filename, vector<complex<double>>& k_s, vector<complex<double>>& offsets){
    fstream file;
    file.open(filename, ios::in);

    string s; getline(file, s); //Skip first line

    while(getline(file, s)){

	if(s.empty())
	    continue;
	//Read off soliton parameters from line
	stringstream ss(s);
	double k_re, o_re, k_im, o_im;
	if(!(ss >> k_re >> k_im >> o_re >> o_im))
	    continue;

	//Create complex number from real and imag part
	complex<double> k(k_re, k_im);
	complex<double> o(o_re, o_im);

	//Store k and conj(k) in vector
	k_s.emplace_back(k);
	k_s.emplace_back(conj(k));

	//Same for o
	offsets.emplace_back(o);
	offsets.emplace_back(conj(o));
    }
}

//Evolve soliton offsets for a fixed time
vector<complex<double>> evolveOffsets(const double t, const vector<complex<double>>& k_s, const vector<complex<double>>& offsets){
    if(k_s.size()!=offsets.size()) throw std::runtime_error("Nonequal k_s and offsets");

    //Initialize new_offset vector
    std::vector<complex<double>> new_offsets;
    new_offsets.reserve(k_s.size());

    int s = k_s.size();
    for(int i = 0; i < s; i++){
	//New offset at time t as function of intial one
	complex<double> new_o = offsets[i] - 12. * k_s[i]*k_s[i]*t;
	new_offsets.emplace_back(new_o);
    }

    return new_offsets;
}

constexpr unsigned long max_points=2*N;
int addPointsToData(const std::vector<Point>& cpoints, std::vector<double>& data){
    //We add the critical points to data.
    //If there are less than max_points, (To be adapted maybe later on for larger N), we add a padding
    //If more than max_points, just skip.
    if(cpoints.size() > max_points){
	cout << "Too many points for writing to numpy" << endl;
	return 0;
    }
    int emplaced = 0;

    for(const Point& p : cpoints){
	data.emplace_back(p.x); data.emplace_back(p.y);
	emplaced++;
    }
    //Now add dummy points
    for(int i = 0; (i + emplaced) < (int)(max_points); i++){
	data.emplace_back(std::numeric_limits<double>::quiet_NaN());
	data.emplace_back(std::numeric_limits<double>::quiet_NaN());
    }
    return 1;
}


int main(){

    //Input/output files
    const string paramfile = "global.params";
    const string solitonfile = "solitons.dat";
    const string output_dir = "Output/";

    //Generate random solitons?
    const bool use_random = false;

    //Read config
    auto config = Config(paramfile);


    auto k_s = vector<complex<double>>{};
    auto offsets = vector<complex<double>>{};


    cout << "Take a moment to check parameters please. " << endl;
    config.Print();

    //Get soliton parameters
    if(use_random){
	//---PARAMETERS FOR RANDOM SOLITON GEN---
	//---COULD BE IN A SEPERATE CONFIG FILE--
	double x_min = -1.; double x_max = 1.;
	double y_min = -1.; double y_max = 1.;
	double max_mod_k = 1.;
	//---------------------------------------

	//Initialize random number generator
	std::random_device dev;
	myRNG rng(dev());

	//Generate parameters
	k_s = randomKs(rng, max_mod_k);
	offsets = randomOffsets(rng, x_min, x_max, y_min, y_max);
    }
    else{
	//Read soliton parameters from file
	readParameters(solitonfile, k_s, offsets);
    }

    //----PRINT OUT INITIAL SOLITON PARAMS-------
    for(const complex<double> & k : k_s)
	cout << k << ", ";
    cout << endl;
    char dummy;
    cin >> dummy;
    //-------------------------------------------

    //Check that number of poles is 
    //same as number in global.hpp
    int pole_number = k_s.size();
    assert(pole_number == N);

    //Define t params
    double current_time = config.t_start;
    double final_time = config.t_end;
    int steps = config.t_points;
    double t_step = (final_time - current_time)/steps;

    //vector for writing critical pts
    std::vector<double> data;
    data.reserve(2*max_points*steps + 1); // A priori two coordinates, and max_points critical pts per step
    unsigned long written = 0;

    //Loop over t values
    while(current_time < config.t_end){

	cout << "Time = " << current_time << endl;

	//Evolve intial offsets for current time
	std::cout << "Evolving offsets" << endl;
	auto t_offsets = evolveOffsets(current_time, k_s, offsets);
	
	//Initialize poles
	cout << "Fetching initial condition" << endl;
	auto [p, v] = initialConditions(k_s, t_offsets, config.y_start);
	PoleState poles(config.y_start, p, v);
	
	//Vector to save at each y
	int total_pt_estimate = config.y_points;
	std::vector<std::unique_ptr<SavedState>> states;
	states.reserve(1.5*total_pt_estimate);

	//Initialize y loop parameters
	double y_step = config.cms_step;
	double y_write_step = (config.y_end - config.y_start)/config.y_points;
	double current_y = config.y_start;

	double y_write_next = config.y_start;
	
	cout << "Starting y loop" << endl;
	Eigen::VectorXd x_vals = Eigen::VectorXd::LinSpaced(config.x_points, config.x_start, config.x_end);
	while(current_y < config.y_end){
	    if(current_y < (y_write_next + 1e-6)){
		poles.insertInto(states);
		y_write_next += y_write_step;
	    }
	    poles.evolveRK4(y_step);
	    current_y = poles.y;
	}
	cout << "Finished CMS evolution." << endl;

	//Finished evolution, now find critical points
	std::vector<Point> critical_points;
	critical_points.reserve(6);

	Eigen::VectorXd prev_phi_x, prev_phi_y;
	cout << "Now looking for critical pts" << endl; 
	if (!states.empty()) {
	    prev_phi_x = phiX(*states[0], config.x_start, config.x_end, config.x_points);
	    prev_phi_y = phiY(*states[0], config.x_start, config.x_end, config.x_points);
	}

	for (std::size_t i = 0; i + 1 < states.size(); ++i) {
	    SavedState& current = *states[i];
	    SavedState& next    = *states[i + 1];

	    // reuse prev_phi_x/prev_phi_y for current
	    Eigen::VectorXd phi_x_next = phiX(next, config.x_start, config.x_end, config.x_points);
	    Eigen::VectorXd phi_y_next = phiY(next, config.x_start, config.x_end, config.x_points);

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
	current_time += t_step;
    }

    //Writing to .npy
    string write_path = output_dir + "cpoints.npy";
    std::vector<unsigned long> shape{written, max_points, 2};
    const npy::npy_data_ptr<double> data_ptr{data.data(), shape, false};
    write_npy(write_path, data_ptr);
}
