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
double Ps, NT;
int constrain;

tuple<int, double *> kin1cbm(double kappa, double X0[2], double Ns);

//driver function
int main() {

	double kmin, kmax, dk, Ns;

	srand(time(0));

	cout << "Enter values for Ps and NT :" << endl;
	cin >> Ps >> NT;

	cout << "Enter kappa range (min, max, dk) :" << endl;
	cin >> kmin >> kmax >> dk;

	while (dk==0)
	{
		cout << "Please enter a non-zero dk :" << endl;
		cin >> dk;
	}

	int numKappas = (kmax-kmin)/dk+1;
	
	if (numKappas == 1)
	{

		cout << "Enter asymptotic sample size : " << endl;
		cin >> Ns;

		double X0[2], X00[2], nSamples;
		int sleepCtr;

		cout << "Enter sample size :" << endl;
		cin >> nSamples;
		
		cout << "Use CBM constraints?" << endl;
		cin >> constrain;

		std::string filename = "bifurc_kappa=" + std::to_string(kmin) + "_Ps=" + std::to_string(Ps) + ".dat";  
		ofstream datafile;
		datafile.open(filename);

		auto start = high_resolution_clock::now();


		for (int n = 0; n < (int) nSamples; n++)
		{
			sleepCtr = 0;
			X00[0] = (double) rand()/RAND_MAX;
			X00[1] = (double) rand()/RAND_MAX;

			while (X00[0]+X00[1]>1)
			{
				X00[0] = (double) rand()/RAND_MAX;
				X00[1] = (double) rand()/RAND_MAX;
			}

			X0[0] = X00[0];
			X0[1] = X00[1];

			auto [t,X] = kin1cbm(kmin,X0,Ns);

			if (t == NT)
			{
				datafile << X00[0] << " " << X00[1]<< " :";

				for (int i = 0; i < Ns; i++)
				{
					datafile << " " << X[i];
				}

				datafile << "\n";
			}

			while (t != NT)
			{
				X00[0] = (double) rand()/RAND_MAX;
				X00[1] = (double) rand()/RAND_MAX;

				while (X00[0]+X00[1]>1)
				{
					X00[0] = (double) rand()/RAND_MAX;
					X00[1] = (double) rand()/RAND_MAX;
				}

				X0[0] = X00[0];
				X0[1] = X00[1];

				auto [t,X] = kin1cbm(kmin,X0,Ns);

				sleepCtr += 1;

				if (t == NT)
				{
					datafile << X00[0] << " " << X00[1]<< " :";

					for (int i = 0; i < Ns; i++)
					{
						datafile << " " << X[i];
					}

					datafile << "\n";
					break;
				}
				else if (sleepCtr==1e6)
				{
					cout << "No valid ICs found after " << std::to_string(sleepCtr) << " attempts." << endl;
					break;
				}
			}
		}

		datafile.close();
		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<milliseconds>(stop - start);

		cout << "Elapsed time : " << duration.count() << " ms" << endl;
	}
	else
	{
		double kappas[numKappas];

		kappas[0] = kmin;
		for (int j = 1; j<numKappas; j++)
		{
			kappas[j] = kappas[j-1]+dk;
		}

		cout << "Enter asymptotic sample size : " << endl;
		cin >> Ns;

		double X0[2], X00[2];
		int singleX0, sleepCtr;

		cout << "Use single X0 ?" << endl;
		cin >> singleX0;

		if (singleX0 == 1)
		{
			cout << "Enter the initial condition X0 :" << endl;
			cin >> X0[0] >> X0[1];

			cout << "Use CBM constraints?" << endl;
			cin >> constrain;

			std::string filename = "bifurc_X0=" + std::to_string(X0[0]) + "_" + std::to_string(X0[1])+ "_Ps=" + std::to_string(Ps) + ".dat";  
			ofstream datafile;
			datafile.open(filename);

			auto start = high_resolution_clock::now();

			for (int k = 0; k < numKappas; k++)
			{
				auto [t,X] = kin1cbm(kappas[k],X0,Ns);

				datafile << kappas[k] << " :";

				for (int i = 0; i < Ns; i++)
				{
					datafile << " " << X[i];
				}

				datafile << "\n";
			}
			datafile.close();
			auto stop = high_resolution_clock::now();
			auto duration = duration_cast<milliseconds>(stop - start);

			cout << "Elapsed time : " << duration.count() << " ms" << endl;
		}
		else if (singleX0 == 0)
		{
			cout << "Use CBM constraints?" << endl;
			cin >> constrain;

			std::string filename = "bifurc_multiX0.dat";  
			ofstream datafile;
			datafile.open(filename);

			auto start = high_resolution_clock::now();

			for (int k = 0; k <= numKappas; k++)
			{
				sleepCtr = 0;
				X00[0] = (double) rand()/RAND_MAX;
				X00[1] = (double) rand()/RAND_MAX;

				while (X00[0]+X00[1]>1)
				{
					X00[0] = (double) rand()/RAND_MAX;
					X00[1] = (double) rand()/RAND_MAX;
				}

				X0[0] = X00[0];
				X0[1] = X00[1];

				auto [t,X] = kin1cbm(kappas[k],X0,Ns);

				if (t == NT)
				{
					datafile << kappas[k] << " " << X00[0] << " " << X00[1]<< " :";

					for (int i = 0; i < Ns; i++)
					{
						datafile << " " << X[i];
					}

					datafile << "\n";
				}

				while (t != NT)
				{
					X00[0] = (double) rand()/RAND_MAX;
					X00[1] = (double) rand()/RAND_MAX;

					while (X00[0]+X00[1]>1)
					{
						X00[0] = (double) rand()/RAND_MAX;
						X00[1] = (double) rand()/RAND_MAX;
					}

					X0[0] = X00[0];
					X0[1] = X00[1];

					auto [t,X] = kin1cbm(kappas[k],X0,Ns);

					sleepCtr += 1;

					if (t == NT)
					{
						datafile << kappas[k] << " " << X00[0] << " " << X00[1]<< " :";

						for (int i = 0; i < Ns; i++)
						{
							datafile << " " << X[i];
						}

						datafile << "\n";
						break;
					}
					else if (sleepCtr==1e6)
					{
						cout << "No valid ICs found after " << std::to_string(sleepCtr) << " attempts." << endl;
						break;
					}
				}
			}

			datafile.close();
			auto stop = high_resolution_clock::now();
			auto duration = duration_cast<milliseconds>(stop - start);

			cout << "Elapsed time : " << duration.count() << " ms" << endl;
		}
	}

	return 0;
}


//CBM iteration
tuple<int, double *> kin1cbm(double kappa, double X0[2], double Ns){

	double * X = new double[(int) Ns];
	double * Z;
	double temp, sumx;
	
	int t, m = 0;

	Z = X0;

	for (t = 1; t<NT; t++)
	{
		temp = (1-Z[0]-Z[1])*(Ps+kappa*(1-Ps)*Z[0]);
		sumx = temp+Z[0];

		if (temp >= 0 && temp <= 1 && sumx <= 1 && constrain==1)
		{
			Z[0] = temp;
			Z[1] = sumx-temp;

			if ((NT-t) <= Ns)
			{
				X[m] = Z[0];
				m += 1;
			}
		}		
		else if (isinf(temp)==0 && isinf(sumx)==0 && constrain==0)
		{
			Z[0] = temp;
			Z[1] = sumx-temp;

			if ((NT-t) <= Ns)
			{
				X[m] = Z[0];
				m += 1;
			}
		}
		else	
		{
			break;
		}
	}

	return {t,X};
}