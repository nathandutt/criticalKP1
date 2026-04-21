#include <array>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
using state = std::array<double, 4>;
using param = std::array<double, 2>;

//Calculate 1/Z^3 force
// s is a vector {re(z), im(z), d/dt re(z), d/dt im(z)}
// and d/dt s = f(s) 
state f(const state& s){

    //|z|^2
    double r2 = s[0]*s[0]+s[1]*s[1];

    return state{s[2], s[3],
	    4.*(s[0]*s[0]*s[0] - 3*s[0]*s[1]*s[1])/(r2*r2*r2),
	    4.*(s[1]*s[1]*s[1]-3*s[1]*s[0]*s[0])/(r2*r2*r2)};
}

//----------Arithmetic operations on state arrays-----------------------

state add(const state& s1, const state& s2){
    return state{s1[0]+s2[0],s1[1]+s2[1],s1[2]+s2[2],s1[3]+s2[3]};
}
state operator+(const state& s1, const state& s2){
    return state{s1[0]+s2[0],s1[1]+s2[1],s1[2]+s2[2],s1[3]+s2[3]};
}

state operator*(const double l, const state& s){
    return state{l*s[0], l*s[1], l*s[2], l*s[3]};
}
double max(double a, double b){
    return (a>b)?a:b;
}

//------------------------------------------------------------------------

//Main RK4 evolution block, with adaptative timestep
double evolveRK4(state& s, const double timestep){
    state k1 = f(s);

    //set adaptive timestep
    double max_f = k1[2]*k1[2] + k1[3]*k1[3];
    double astep = (max_f<1.) ? timestep : timestep / std::pow(max_f, 0.3);

    //Calculate successive terms for RK4
    state k2 = f(s + 0.5*astep*k1);
    state k3 = f(s + 0.5*astep*k2);
    state k4 = f(s + astep*k3);

    //Update s
    s = s + (astep/6.) * (k1 + 2.*k2 + 2.*k3 + k4);

    return astep;
}


//initial state defined by an offset, that will depend on real time, and delta k
state initialState(const param& k1, const param& k2, const double yi, const double t){
    double toffs_re = -12.*(k1[0]*k1[0]-k1[1]*k1[1] - (k2[0]*k2[0]-k2[1]*k2[1]));
    double toffs_im = -12.*(2.*k1[0]*k1[1] - 2.*k2[0]*k2[1]);
    return state{2.*(k1[0] - k2[0])*yi + toffs_re*t, 2.*(k1[1] - k2[1])*yi + toffs_im*t,
		 2.*(k1[0] - k2[0]), 2.*(k1[1] - k2[1])};
}
state averageState(const param& k1, const param& k2, const double y, const double t){
    double toffs_re = -12.*(k1[0]*k1[0]-k1[1]*k1[1] + (k2[0]*k2[0]-k2[1]*k2[1]));
    double toffs_im = -12.*(2.*k1[0]*k1[1] + 2.*k2[0]*k2[1]);
    return state{2.*(k1[0] + k2[0])*y + toffs_re*t, 2.*(k1[1] + k2[1])*y + toffs_im*t,
		 2.*(k1[0] + k2[0]), 2.*(k1[1] + k2[1])};
}

int main(){

    //Define simulation parameters
    double ti = -2;
    double tf = 2.;
    double tcurr = ti;
    double tstep = 0.1;

    double yi = -50;
    double yf = 50.;
    double ystep = 0.01;

    //Define particle parameters
    param k2{0., 0.6};
    param k1{0., 0.3};
    param offs{0., 0.};

    //Open output files
    fstream file1;
    fstream file2;
    file1.open("output1.res", ios::out);
    file2.open("output2.res", ios::out);

    //Run integration
    while(tcurr < tf){
	std::cout << "T = " << tcurr << endl;
	tcurr += tstep;
	double ycurr = yi;
	state s = initialState(k1, k2, yi, tcurr);
	double ywrite = yi;
	while (ywrite < yf) {
	    if (ywrite > yf) {
		ywrite = yf;
	    }

	    while (ycurr < ywrite) {
		double step = evolveRK4(s, ystep);
		ycurr += step;
		
	    }
	    state z_avg = averageState(k1, k2, ycurr, tcurr);
	    state s1 = 0.5*add(z_avg, s);
	    state s2 = 0.5*add(z_avg, -1.*s);
	    file1 << s1[0] << " " << s1[1] << " ";
	    file2 << s2[0] << " " << s2[1] << " ";
	    ywrite += ystep;
	}
	file1 << endl;
	file2 << endl;
    }
    file1.close();
    file2.close();
}
