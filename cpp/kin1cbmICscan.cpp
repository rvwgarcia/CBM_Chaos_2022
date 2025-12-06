/*
	[UPDATE FROM cxmodelChaosICs.cpp: Designed to produce IC admissibility data
	for ML algorithms. 10k samples (Kappa=3.5, Ps=0, NT=1e6) takes about 6 min.]

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
#include <ctime>
#include <tuple>
#include <math.h>
#include <fstream>

#include <omp.h>
#include <unistd.h>

using namespace std;
using namespace std::chrono;

const int taur = 2;
double NT;

void icGen(double *X0);
tuple<int, double, double> kin1cbm(double kappa, double Ps, double X0[2]);

//driver function
int main()
{
	double nSamples;
	int FixedKappa;
	srand(time(0));//^ (1+omp_get_thread_num())

	cout << "Fixed kappa ?" << endl;
	cin >> FixedKappa;

	//cout << "Enter values for Ps and NT :" << endl;
	//cin >> Ps >> NT;
	NT = 1e6;

	if (FixedKappa==0)
	{
		cout << "Enter number of samples :" << endl;
		cin >> nSamples;

		auto start = high_resolution_clock::now();
		std::time_t start_time = high_resolution_clock::to_time_t(start);
		string filename = std::ctime(&start_time);

		ofstream datafile;
		datafile.open("ScanBCs_"+filename+".dat");

		#pragma omp parallel for
		for (int i = 0; i < (int) nSamples; i++)
		{	
			double X0[2];

			usleep(5000 * omp_get_thread_num());
			
			double kappa = (3.6760-1)*(double) rand()/RAND_MAX+1;
			double Ps = (double) rand()/RAND_MAX;

			icGen(X0);

			auto [t,mle,cumulant] = kin1cbm(kappa,Ps,X0);

			#pragma omp critical
			{
				datafile << X0[0] << " " << X0[1] << " " << kappa << " " << Ps << " " << mle << " " << cumulant << " " << t << "\n";
			}
		}

		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<milliseconds>(stop - start);
		cout << "Elapsed time : " << duration.count() << " ms" << endl;

		datafile.close();
	}
	else if (FixedKappa==1)
	{
		double kappa, Ps=0;

		cout << "Enter kappa value :" << endl;
		cin >> kappa;

		cout << "Enter number of samples :" << endl;
		cin >> nSamples;

		ofstream datafile;
		std::string filename = "ScanICs_Fixed-Kappa" + std::to_string(kappa) + ".dat";
		datafile.open(filename);

		auto start = high_resolution_clock::now();

		#pragma omp parallel for
		for (int i = 0; i < (int) nSamples; i++)
		{	
			//int sleepCtr = 0;
			double X0[2];

			usleep(5000 * omp_get_thread_num());

			//X0[2] = icGen();
			icGen(X0);
			//cout << X0[0] << " " << X0[1] << endl;
			auto [t,mle,cumulant] = kin1cbm(kappa,Ps,X0);
			//cout << X0[0] << " " << X0[1] << endl;
			//cout << X[0] << " " << X[1] << endl;
			#pragma omp critical
			{
				datafile << X0[0] << " " << X0[1] << " " << kappa << " " << Ps << " " << mle << " " << cumulant << " " << t << "\n";
			}
		}

		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<milliseconds>(stop - start);
		cout << "Elapsed time : " << duration.count() << " ms" << endl;

		datafile.close();
	}

	return 0;
}

void icGen(double X0[2])
{
	X0[0] = (double) rand()/RAND_MAX;
	X0[1] = (double) rand()/RAND_MAX;

	while (X0[0]+X0[1]>1)
	{
		X0[0] = (double) rand()/RAND_MAX;
		X0[1] = (double) rand()/RAND_MAX;	//[omp_get_thread_num()]
	}
}

//CBM iteration and Lyapunov calculation
tuple<int, double, double> kin1cbm(double kappa, double Ps, double X0[2])
{
	double X[2];
	double logdX[(int) NT];

	double A[2][2], v[2], dX[2];
	double ProdDF[2][2] = {{1,0},{0,1}};

	double c = kappa*(1-Ps), temp, temp2, sum, sumx, norm0, normn, cumulant = 0, mle = 0, sum2 = 0, sum4 = 0;

	int t, m = 0, breaker = 0;

	X[0] = X0[0];
	X[1] = X0[1];

	for (t = 1; t<NT; t++)
	{
		//This part prepares data to calculate the Lyapunov exponent and cumulant

		if (breaker==0) //((NT-t) <= 1e4 && breaker==0)
		{
			dX[0] = (double) rand()/(1e6*RAND_MAX);
			dX[1] = (double) rand()/(1e6*RAND_MAX);

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
				breaker = 0;
			}
			else
			{
				ProdDF[0][0] = A[0][0];
				ProdDF[1][0] = A[1][0];
				ProdDF[0][1] = A[0][1];
				ProdDF[1][1] = A[1][1];

				logdX[m] = log(normn/norm0)/(m+1);

				sum += logdX[m];
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
			break;
		}
	}

	if (m > 0)
	{
		mle = logdX[m-1];
		//mle = sum / m;
		for (int n = 0; n < m; n++)
		{
			temp2 = logdX[n]-mle;
			sum2 += pow(temp2,2);
			sum4 += pow(temp2,4);
		}
		
		cumulant = sum4/m-3*pow(sum2/m,2);
	}

	return {t, mle, cumulant};
}