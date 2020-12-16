#define _USE_MATH_DEFINES
#include <iostream>
#include <random>
#include "parameters.h"
using namespace std;

random_device seed;
mt19937 generator(seed());
uniform_real_distribution<double> uniform(0, 1);
exponential_distribution<double> exponential(1);
normal_distribution<double> normal(0, 1);



void EnergyProfile(double ***chCoef)
{
	for(int resource = 0; resource < 2; resource++){
		for (int nuser = 0; nuser < NUM_USER * SCMA_SOURCE / SCMA_USER_SOURCE; nuser++)
		{
			double theta = 2 * M_PI * uniform(generator);		// phase shift
			chCoef[resource][nuser][4] = theta;
			chCoef[resource][nuser][3] = exponential(generator);			// energy
			chCoef[resource][nuser][2] = sqrt(chCoef[resource][nuser][3]);			// amplitude
			chCoef[resource][nuser][1] = chCoef[resource][nuser][2] * sin(theta);	// imaginary pary
			chCoef[resource][nuser][0] = chCoef[resource][nuser][2] * cos(theta);	// real part


		//	cout << theta << " " << chCoef[nuser][2] << endl;
		}
	}
}

void MultipleAccessChannel(double stdDev, double ***chCoef, double ***chCoef2, double **tx, double **rx, int index, int scma_matrix[SCMA_SOURCE][SCMA_SOURCE * NUM_USER / SCMA_USER_SOURCE], int coef_idx[SCMA_SOURCE][SCMA_SOURCE * NUM_USER / SCMA_USER_SOURCE])
{
	int nresource;
	for (int i = 0; i < BLOCK_LEN; ++i)
	{
		rx[i][0] = stdDev * normal(generator); // real
		rx[i][1] = stdDev * normal(generator); // imaginary
		for (int nuser = 0; nuser < NUM_USER * SCMA_SOURCE / SCMA_USER_SOURCE; ++nuser)
		{

			if (scma_matrix[index][nuser] == 1) {
				nresource = coef_idx[index][nuser]; 
				rx[i][0] += tx[nuser][2 * i + 0] * chCoef[nresource][nuser][0];
				rx[i][0] -= tx[nuser][2 * i + 1] * chCoef2[nresource][nuser][1];
				rx[i][1] += tx[nuser][2 * i + 0] * chCoef[nresource][nuser][1];
				rx[i][1] += tx[nuser][2 * i + 1] * chCoef2[nresource][nuser][0];
			}
		}
	}
}