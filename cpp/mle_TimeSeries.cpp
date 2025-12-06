/*
	This code may take a while to run.

	LyapStats produces a statistical average of characteristic Lyapunov exponents
	(MLEs) by sampling different ICs at each kappa value.

	With LyapStats turned off, the code calculates the MLE over a range of kappa
	using a single IC. This calculation is performed nSamples number of times with
	nSamples different ICs.

	If nSamples=1, the code will calculate the MLE over the range of kappa values
	given.

	\rvwg2022
*/

#include <iostream>
#include <stdlib.h>
#include <chrono>
#include <tuple>
#include <math.h>
#include <fstream>

#include <omp.h>
#include <unistd.h>

using namespace std;
using namespace std::chrono;

const int taur = 2;
double Ps, NT, kappa;

//driver function
int main() {
	srand(time(0));//^ (1+omp_get_thread_num())

	//cout << "Enter values for kappa, Ps, and NT :" << endl;
	//cin >> kappa >> Ps >> NT;
	kappa = 3.65;
	Ps = 0;
	NT = 1E6;

	double X0[2];
	//cout << "Enter the initial condition X0 :" << endl;
	//cin >> X0[0] >> X0[1];
	X0[0] = 0.3100;
	X0[1] = 0.1211;

	ofstream datafile;
	auto start = high_resolution_clock::now();

 
	std::string filename = "mle_Series_X0=" + std::to_string(X0[0]) + "_" + std::to_string(X0[1])+ "_k=" + std::to_string(kappa) + "_Ps=" + std::to_string(Ps) + ".dat";
	datafile.open(filename);

	double X[2];
	double logdX;

	double A[2][2], v[2], dX[2];
	double ProdDF[2][2] = {{1,0},{0,1}};

	double c = kappa*(1-Ps), temp, temp2, sumx, norm0, normn;

	int t, m = 0, breaker = 0;
	double divConstant = RAND_MAX;

	X[0] = X0[0];
	X[1] = X0[1];

	for (t = 1; t<NT; t++)
	{
		//This part prepares data to calculate the Lyapunov exponent and cumulant

		if  (breaker==0)//(NT-t) <= 1e4 && 
		{
			dX[0] = (double) rand()/divConstant;
			dX[1] = (double) rand()/divConstant;

			norm0 = sqrt(dX[0] * dX[0] + dX[1] * dX[1]);
			// A = DF(t)*ProdDF;
			A[0][0] = ((1-X[0]-X[1])*c-(Ps+c*X[0])) * ProdDF[0][0]-(Ps+c*X[0]) * ProdDF[1][0];
			A[0][1] = ((1-X[0]-X[1])*c-(Ps+c*X[0])) * ProdDF[0][1]-(Ps+c*X[0]) * ProdDF[1][1];
			A[1][0] = ProdDF[0][0];
			A[1][1] = ProdDF[0][1];

			// v = A*dX0;
			v[0] = A[0][0] * dX[0] + A[0][1] * dX[1];
			v[1] = A[1][0] * dX[0] + A[1][1] * dX[1];
			normn = sqrt(v[0] * v[0] + v[1] * v[1]);

			if (normn==0 || isinf(normn)==1)
			{
				cout << "Vector norm singularity (norm = " << normn << ") encountered after " << t << " timesteps." << endl;
				breaker = 0;
			}
			else
			{
				ProdDF[0][0] = A[0][0];
				ProdDF[1][0] = A[1][0];
				ProdDF[0][1] = A[0][1];
				ProdDF[1][1] = A[1][1];

				logdX = log(normn/norm0)/(m+1);
				cout << logdX << endl;
				datafile << t << " " << logdX << "\n";
				m += 1;
			}

			
		}

		temp = (1-X[0]-X[1])*(Ps+kappa*(1-Ps)*X[0]);
		sumx = temp+X[0];

		if (temp >= 0 && temp <= 1 && sumx <= 1)
		{
			X[0] = temp;
			X[1] = sumx-temp;
		}
		else
		{
			cout << "Inadmissible initial condition encountered!" << endl;
			break;
		}
	}

		

	datafile.close();

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);

	cout << "Elapsed time : " << duration.count() << " ms" << endl;

	return 0;
}