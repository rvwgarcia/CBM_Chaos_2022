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



double * icGen();

int main(){
	//int nThreads = omp_get_max_threads();
	int nSamples = 8;
	int m = 0;

	//double X0[nThreads][2];
	double * X0;

	#pragma omp parallel for
	for (int i = 0; i < nSamples; i++)
	{
		usleep(500 * omp_get_thread_num());
		m = omp_get_thread_num();
		X0 = icGen();

		cout << m << " " << X0[0] << endl;
		//X0[m][0] = temp[0];
		//X0[m][1] = temp[1];
		//cout << m << " " << temp[0] << endl;
	}

	return 0;
}

double * icGen(){
	double * X0 = new double[2];

	X0[0] = (double) rand()/RAND_MAX;
	X0[1] = (double) rand()/RAND_MAX;

	while (X0[0]+X0[1]>1)
	{
		X0[0] = (double) rand()/RAND_MAX;
		X0[1] = (double) rand()/RAND_MAX;	//[omp_get_thread_num()]
	}
	return X0;
}