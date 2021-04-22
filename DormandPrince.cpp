#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <array>

using namespace std;

//integrator parameters
struct integparam {             
	double h = 0.00005;         // initial step, in seconds
	double tolerance = 1e-3;    // integrator tolerance
	double t0 = 0;              // start time
	double t_end = 3;           // end time
	double y0 = 1;
};

std::array<double, 7> c7 = {
	0.e0,      1.e0 / 5.e0,       3.e0 / 10.e0,        4.e0 / 5.e0,      8.e0 / 9.e0,           1e0,       1.e0 };

static const std::array<double, 49> a49 = {
	0.e0,      0.e0,      0.e0,      0.e0,      0.e0,      0.e0,      0.e0,
	1.e0 / 5.e0,      0.e0,      0.e0,      0.e0,      0.e0,      0.e0,      0.e0,
	3.e0 / 40.e0,	  9.e0 / 40.e0,	     0.e0,       0.e0,       0.e0,       0.e0,       0.e0,
	44.e0 / 45.e0,      -56.e0 / 15.e0,      32.e0 / 9.e0,       0.e0,        0.e0,        0.e0,        0.e0,
	19372.e0 / 6561.e0,      -25360.e0 / 2187.e0,      64448.e0 / 6561.e0,       -212.e0 / 729.e0,        0.e0,        0.e0,        0.e0,
	9017.e0 / 3168.e0,      -355.e0 / 33.e0,      46732.e0 / 5247.e0,       49.e0 / 176.e0,        -5103.e0 / 18656.e0,        0.e0,        0.e0,
	35.e0 / 384.e0,       0.e0,      500.e0 / 1113.e0,       125.e0 / 192.e0,        -2187.e0 / 6784.e0,        11.e0 / 84.e0,        0.e0 };

static const std::array<double, 7> b27 = {
	35.e0 / 384.e0,       0.e0,      500.e0 / 1113.e0,       125.e0 / 192.e0,        -2187.e0 / 6784.e0,        11.e0 / 84.e0,        0.e0 };

static const std::array<double, 7> b17 = {
	5179.e0 / 57600,      0.e0,       7571.e0 / 16695.e0,        393.e0 / 640.e0,       -920997.e0 / 339200.e0,           187.e0 / 2100.e0,       1.e0 / 40.e0 };

double func(double t,double y) {
	return y*cos(3*t); //equation right side
}

typedef std::array<double, 2> OutArray;

void DormandPrince(double t0, double y0, double h, OutArray* out)
{
	std::array<double, 7> k;
	//this is to compute the values of k
	k[0] = func(t0, y0);
	k[1] = func(t0 + c7[1] * h, y0 + h * a49[7] * k[0]);
	k[2] = func(t0 + c7[2] * h, y0 + h * (a49[14] * k[0] + a49[15] * k[1]));
	k[3] = func(t0 + c7[3] * h, y0 + h * (a49[21] * k[0] + a49[22] * k[1] + a49[23] * k[2]));
	k[4] = func(t0 + c7[4] * h, y0 + h * (a49[28] * k[0] + a49[29] * k[1] + a49[30] * k[2] + a49[31] * k[3]));
	k[5] = func(t0 + c7[5] * h, y0 + h * (a49[35] * k[0] + a49[36] * k[1] + a49[37] * k[2] + a49[38] * k[3] + a49[39] * k[4]));
	k[6] = func(t0 + c7[6] * h, y0 + h * (a49[42] * k[0] + a49[43] * k[1] + a49[44] * k[2] + a49[45] * k[3] + a49[46] * k[4] + a49[47] * k[5]));
	(*out)[0] = y0 + h * (b17[0] * k[0] + b17[1] * k[1] + b17[2] * k[2] + b17[3] * k[3] + b17[4] * k[4] + b17[5] * k[5] + b17[6] * k[6]);
	(*out)[1] = y0 + h * (b27[0] * k[0] + b27[1] * k[1] + b27[2] * k[2] + b27[3] * k[3] + b27[4] * k[4] + b27[5] * k[5] + b27[6] * k[6]);
}

int main() {
	integparam integr;
	OutArray out;

	double est;
	double t = integr.t0;
	double h = integr.h;
	double y0 = integr.y0;
	std::cout.setf(std::ios::fixed);
	std::ofstream  resultFile("result.csv");
	
	while (t < integr.t_end)
	{
		DormandPrince(t, y0, h, &out);
		est = abs(out[0] - out[1]);
		if (est > integr.tolerance) {
			h = h * pow((integr.tolerance / est), 1.0 / 5.0);
			cout << h << endl;
		}
		else {
			//cout << "1" << endl;
			resultFile << t << " " << out[1] << endl;
			y0 = out[1];
			t = t + h;
		}
	}

	h = integr.t_end - t + h;
	DormandPrince(t, y0, h, &out);
	resultFile << t << " " << out[1] << endl;
	resultFile.close();

	return 0;
}
