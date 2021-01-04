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
		chCoef[nuser][3] = exponential(generator);			// energy
		chCoef[nuser][2] = sqrt(chCoef[nuser][3]);			// amplitude
		chCoef[nuser][1] = chCoef[nuser][2] * sin(theta);	// imaginary pary
		chCoef[nuser][0] = chCoef[nuser][2] * cos(theta);	// real part
		if((nuser % 2) == 1){
			chCoef[nuser][3] = chCoef[nuser - 1][3]	;		// energy
			chCoef[nuser][2] = sqrt(chCoef[nuser][3]);			// amplitude
			chCoef[nuser][1] = chCoef[nuser - 1][0];	// imaginary pary
			chCoef[nuser][0] = -chCoef[nuser - 1][1];	// real part
		}
	}
}

void MultipleAccessChannel(double stdDev, double **chCoef, double **tx, double **rx,double **pilot, int *known_drift)
{
	int effLen;
	if(COLLISION)
		effLen = SYNCHRONOUS ? BLOCK_LEN : BLOCK_LEN * UP_RATE * SPREAD_LEN;
	else
		effLen = SYNCHRONOUS ? known_drift[NUM_USER-1]+BLOCK_LEN : (known_drift[NUM_USER - 1] + BLOCK_LEN ) * UP_RATE * SPREAD_LEN;

	for (int i = 0; i < effLen; i++)
	{
		rx[i][0] = stdDev * normal(generator); // real
		rx[i][1] = stdDev * normal(generator); // imaginary
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			if (!SYNCHRONOUS)
			{
				rx[i][0] += tx[nuser][i] * chCoef[nuser][0] - (Power_Ratio) / (1 - Power_Ratio) * pilot[nuser][i] * chCoef[nuser][1];
				rx[i][1] += tx[nuser][i] * chCoef[nuser][1] + (Power_Ratio) / (1 - Power_Ratio) * pilot[nuser][i] * chCoef[nuser][0];
			}
			else
			{
				rx[i][0] += tx[nuser][i] * chCoef[nuser][0];
				rx[i][1] += tx[nuser][i] * chCoef[nuser][1];
			}
			tx[nuser][i] = 0;
		}
	}
}