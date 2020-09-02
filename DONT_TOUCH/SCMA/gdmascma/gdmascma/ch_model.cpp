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

void EnergyProfile(double **chCoef)
{
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		double theta = 2 * M_PI * uniform(generator);		// phase shift
		chCoef[nuser][4] = theta;
		chCoef[nuser][3] = exponential(generator);			// energy
		chCoef[nuser][2] = sqrt(chCoef[nuser][3]);			// amplitude
		chCoef[nuser][1] = chCoef[nuser][2] * sin(theta);	// imaginary pary
		chCoef[nuser][0] = chCoef[nuser][2] * cos(theta);	// real part


	//	cout << theta << " " << chCoef[nuser][2] << endl;
	}
}

void MultipleAccessChannel(double stdDev, double **chCoef, double **tx, double **rx)
{
	for (int i = 0; i < BLOCK_LEN; i++)
	{
		rx[i][0] = stdDev * normal(generator); // real
		rx[i][1] = stdDev * normal(generator); // imaginary
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			rx[i][0] += tx[nuser][2 * i + 0] * chCoef[nuser][0];
			rx[i][0] -= tx[nuser][2 * i + 1] * chCoef[nuser][1];
			rx[i][1] += tx[nuser][2 * i + 0] * chCoef[nuser][1];
			rx[i][1] += tx[nuser][2 * i + 1] * chCoef[nuser][0];
		}
	}
}