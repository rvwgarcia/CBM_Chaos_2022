#include <iostream>
#include <stdlib.h>
#include <chrono>
#include <math.h>
#include <fstream>

#include <omp.h>
#include <unistd.h>

using namespace std;
using namespace std::chrono;

const int taur = 2;
double NT;

int kin1cbm(double kappa, double Ps, double X[2]);


//driver function
int main() {
	double kmin, kmax, kappa, Ps, nSamples, X0[2], X00[2];
	int maxIters;

	srand(time(0));

	cout << "Enter max iterations :" << endl;
	cin >> NT;

	cout << "Enter kappa range (min, max) :" << endl;
	cin >> kmin >> kmax;

	cout << "Enter number of samples :" << endl;
	cin >> nSamples;

	cout << "Enter the initial condition X0 :" << endl;
	cin >> X00[0] >> X00[1];

	ofstream datafile;
	std::string filename = "Mandlebrot_X0=" + std::to_string(X00[0]) + "_" + std::to_string(X00[1])+ ".dat";
	datafile.open(filename);

	auto start = high_resolution_clock::now();

	#pragma omp for
	for (int n = 0; n<(int) nSamples; n++)
	{
		usleep(5000 * omp_get_thread_num());

		kappa = (double) rand()/RAND_MAX * (kmax-kmin)+kmin;
		Ps = (double) rand()/RAND_MAX;

		X0[0] = X00[0];
		X0[1] = X00[1];

		maxIters = kin1cbm(kappa, Ps, X0);

		datafile << Ps << " " << kappa << " " << maxIters << "\n";
	}

	datafile.close();

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);

	cout << "Elapsed time : " << duration.count() << " ms" << endl;

	return 0;
}


//CBM iteration
int kin1cbm(double kappa, double Ps, double X[2]){
	double temp, sumx;
	int t;

	for (t = 1; t<NT; t++)
	{
		temp = (1-X[0]-X[1])*(Ps+kappa*(1-Ps)*X[0]);
		sumx = temp+X[0];

		if (temp >= 0 && temp <= 1 && sumx <= 1)
		{
			X[0] = temp;
			X[1] = sumx-temp;
		}
		else
		{
			break;
		}
	}

	return t;
}