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

int kin, taur;
double Ps, NT, P1, P2;

double * icGen();
tuple<int, double *, double, double> mfcbm(double kappa, double * X0);

//driver function
int main()
{
	int kin, taur;
	double kmin, kmax, dk, B;
	int nSamples, LyapStats;//, nThreads = omp_get_max_threads();

	srand(time(0));//^ (1+omp_get_thread_num())

	cout << "Enter values for Kin and Taur :" << endl;
	cin >> kin >> taur;

	if (kin==2)
	{
		//cout << "Enter values for B :" << endl;
		//cin >> B;
		B = 0.5;
		P1 = exp(-B) / (exp(-B) + exp(-2*B));
		P2 = exp(-2*B) / (exp(-B) + exp(-2*B));
	}

	cout << "Enter values for Ps and NT :" << endl;
	cin >> Ps >> NT;

	cout << "Enter kappa range (min, max, dk) :" << endl;
	cin >> kmin >> kmax >> dk;

	cout << "Enter number of samples :" << endl;
	cin >> nSamples;

	if (nSamples==1)
	{
		while (dk==0)
		{
			cout << "Please enter a non-zero dk :" << endl;
			cin >> dk;
		}

		double X0[2];
		cout << "Enter the initial condition X0 :" << endl;
		cin >> X0[0] >> X0[1];

		ofstream datafile;
		auto start = high_resolution_clock::now();

		if(kmin == kmax)
		{
			std::string filename = "output_SingleX0=" + std::to_string(X0[0]) + "_" + std::to_string(X0[1])+ "_k=" + std::to_string(kmin) + "_Ps=" + std::to_string(Ps) + ".dat";
			datafile.open(filename);
		}
		else
		{
			std::string filename = "output_SingleX0=" + std::to_string(X0[0]) + "_" + std::to_string(X0[1])+ "_Ps=" + std::to_string(Ps) + ".dat"; 
			datafile.open(filename);
		}

		for (double kappa = kmin; kappa <= kmax; kappa += dk)
		{
			auto [t,X,mle,cumulant] = mfcbm(kappa, X0);

			datafile << kappa << " " << t << " " << mle << " " << cumulant << "\n";
		}

		datafile.close();

		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<milliseconds>(stop - start);

		cout << "Elapsed time : " << duration.count() << " ms" << endl;
	}
	else if (kmin == kmax && nSamples>1)//Parallelized nSamples>1, fixed kappa
	{
		double kappa = kmin;

		ofstream datafile;
		std::string filename = "output_FixedKappa" + std::to_string(kappa) + ".dat";
		datafile.open(filename);

		auto start = high_resolution_clock::now();

		double * X0;

		#pragma omp parallel for
		for (int i = 0; i < nSamples; i++)
		{	
			int sleepCtr = 0;

			usleep(1000 * omp_get_thread_num());

			X0 = icGen();
			auto [t,X,mle,cumulant] = mfcbm(kappa, X0);

			if (t == NT)
			{
				#pragma omp critical
				{
					datafile << X0[0] << " " << X0[1] << " " << mle << " " << cumulant << "\n";
				}
			}
			while (t != NT)
			{
				X0 = icGen();
				auto [t,X,mle,cumulant] = mfcbm(kappa, X0);

				sleepCtr += 1;

				if (t == NT)
				{
					#pragma omp critical
					{
						datafile << X0[0] << " " << X0[1] << " " << mle << " " << cumulant << "\n";
					}
					break;
				}
				else if (sleepCtr==1e6)
				{
					cout << "No valid ICs found after " << std::to_string(sleepCtr) << " attempts." << endl;
					break;
				}
			}
		}

		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<milliseconds>(stop - start);
		cout << "Elapsed time : " << duration.count() << " ms" << endl;

		datafile.close();
	}
	else if(kmin != kmax && nSamples>1)	//Parallelized nSamples>1, kappa sweep
	{
		while (dk==0)
		{
			cout << "Please enter a non-zero dk :" << endl;
			cin >> dk;
		}

		cout << "Perform Lyapunov statistics ?" << endl;
		cin >> LyapStats;

		if (LyapStats == 1)
		{
			double mleAvg, cumulantAvg;
			double * X0;

			auto start = high_resolution_clock::now();
			
			ofstream datafile;
			std::string filename = "LyapStats_Ps=" + std::to_string(Ps) + ".dat";  
			datafile.open(filename);

			int nKappas = (kmax-kmin)/dk+1;

			#pragma omp for
			for (int k = 0; k <= nKappas; k++)
			{
				double kappa = kmin + dk * k;
				
				double mleSum = 0, cumulantSum = 0;

				usleep(1000 * omp_get_thread_num());

				for (int i = 0; i < nSamples; i++)
				{
					int sleepCtr = 0;
					X0 = icGen();

					auto [t,X,mle,cumulant] = mfcbm(kappa, X0);

					if (t == NT)
					{
						mleSum += mle;
						cumulantSum += cumulant;
					}
					while (t != NT)
					{
						X0 = icGen();

						auto [t,X,mle,cumulant] = mfcbm(kappa, X0);

						sleepCtr += 1;

						if (t == NT)
						{
							mleSum += mle;
							cumulantSum += cumulant;
							break;
						}
						else if (sleepCtr==1e6)
						{
							cout << "No valid ICs found after " << std::to_string(sleepCtr) << " attempts." << endl;
							break;
						}
					}
				}

				mleAvg = mleSum / nSamples;
				cumulantAvg = cumulantSum / nSamples;

				datafile << kappa << " " << mleAvg << " " << cumulantAvg << "\n";
			}

			datafile.close();

			auto stop = high_resolution_clock::now();
			auto duration = duration_cast<milliseconds>(stop - start);

			cout << "Elapsed time : " << duration.count() << " ms" << endl;
		}
		else if (LyapStats == 0)
		{
			#pragma omp parallel for
			for (int i = 0; i < nSamples; i++)
			{	
				ofstream datafile;

				auto start = high_resolution_clock::now();
				double * X0;

				X0 = icGen();

				std::string filename = "output_" + std::to_string(X0[0]) + "_" + std::to_string(X0[1]) + ".dat";  

				datafile.open(filename);

				for (double kappa = kmin; kappa <= kmax; kappa += dk)
				{
					auto [t,X,mle,cumulant] = mfcbm(kappa, X0);

					datafile << kappa << " " << t << " " << mle << " " << cumulant << "\n";
				}

				datafile.close();

				auto stop = high_resolution_clock::now();
				auto duration = duration_cast<milliseconds>(stop - start);

				#pragma omp critical
				{
					cout << "Elapsed time : " << duration.count() << " ms" << endl;
				}
			}
		}
	}

	return 0;
}

double * icGen(){
	double * X0 = new double[taur];
	double sumx = 0;

	for (int i=0; i<=taur; i++)
	{
		X0[i] = (double) rand()/RAND_MAX;
		sumx += X0[i];
	}

	while (sumx>1)
	{
		sumx = 0;
		for (int i=0; i<=taur; i++)
		{
			X0[i] = (double) rand()/RAND_MAX;
			sumx += X0[i];
		}
	}

	return X0;
}

//CBM iteration and Lyapunov calculation
tuple<int, double *, double, double> mfcbm(double kappa, double * X0){
	double * X;
	double * logdX = new double[(int) NT];

	double A[taur][taur], v[taur], dX[taur];
	double ProdDF[taur][taur];

	for (int i=0; i<taur; i++)
	{
		for (int j=0; j<taur; j++)
		{
			if (i==j)
				ProdDF[i][j] = 1;
			else
				ProdDF[i][j] = 0;
		}
	}

	double c = kappa*(1-Ps), temp, sumx, sumv, sumAdX, sumProds, norm0, normn, sumdX2;
	double cumulant = 0, mle = 0, sum2 = 0, sum4 = 0;

	int i, j, t, m = 0, breaker = 0;

	X = X0;

	for (t = 1; t<NT; t++)
	{
		//This part prepares data to calculate the Lyapunov exponent and cumulant

		if (breaker==0) //((NT-t) <= 1e4 && breaker==0)
		{
			sumdX2 = 0;
			sumx = 0;
			for (i=0; i<taur; i++)
			{
				dX[i] = (double) rand()/(1e6*RAND_MAX);
				sumdX2 += dX[i] * dX[i];
				sumx += X[i];
			}

			norm0 = sqrt(sumdX2);
			
			// A = DF(t)*ProdDF;
			for (j=0;j<taur;j++)
			{
				sumProds = 0;
				for (i=0;i<taur;i++)
				{
					sumProds += ProdDF[i][j];
				}

				if (kin==1)
				{
					A[0][j] = ((1-sumx)*c-(Ps+c*X[0])) * ProdDF[0][j]-(Ps+c*X[0]) * sumProds;
				}
				else if (kin==2)
				{
					A[0][j] = ((1-sumx)*(c-2*kappa*P1*P2*c*X[0])-(Ps+c*X[0]-kappa*c*P1*P2*X[0]*X[0])) * ProdDF[0][j]-(Ps+c*X[0]-kappa*c*P1*P2*X[0]*X[0]) * sumProds;
				}
			}

			for (i=1;i<taur;i++)
			{
				for (j=0;j<taur;j++)
				{
					A[i][j] = ProdDF[i-1][j];
				}
			}

			// v = A*dX0;
			
			sumv = 0;
			for (i=0; i<taur; i++)
			{
				sumAdX = 0;
				for (j=0; j<taur; j++)
				{
					sumAdX += A[i][j] * dX[i];
				}
				
				v[i] = sumAdX;
				sumv += v[i] * v[i];
			}

			normn = sqrt(sumv);

			if (normn==0 || isinf(normn)==1)
			{
				breaker = 0;
			}
			else
			{
				for (i=0; i<taur; i++)
				{
					for (j=0; j<taur; j++)
					{
						ProdDF[i][j] = A[i][j];
					}
				}

				logdX[m] = log(normn/norm0)/(m+1);
				//sum += logdX[m];

				m += 1;
			}
		}

		sumx = 0;//sum at n
		for (i=0; i<taur; i++)
		{
			sumx += X[i];
		}

		if (kin==1)
		{
			temp = (1-sumx)*(Ps+kappa*(1-Ps)*X[0]);
		}
		else if (kin==2)
		{
			temp = (1-sumx)*(Ps+c*X[0]-kappa*c*P1*P2*X[0]*X[0]);
		}
		
		sumx = temp+sumx-X[taur-1]; //sum at n+1

		if (temp >= 0 && temp <= 1 && sumx <= 1)
		{
			for (i=1;i<taur;i++)
			{
				//X[1] = sumx-temp;
				X[i] = X[i-1];
			}

			X[0] = temp;
		}
		else
		{
			break;
		}
	}

	if (m > 0)
	{
		mle = logdX[m-1];//mle = sum / m;
		
		for (int n = 0; n < m; n++)
		{
			temp = logdX[n]-mle;
			sum2 += pow(temp,2);
			sum4 += pow(temp,4);
		}
		
		cumulant = sum4/m-3*pow(sum2/m,2);
	}

	delete [] logdX;

	return {t, X, mle, cumulant};
}