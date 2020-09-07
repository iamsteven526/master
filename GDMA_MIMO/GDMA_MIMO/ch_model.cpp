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
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		for (int i = 0; i < NUM_TX; i++)
		{
			double theta = 2 * M_PI * uniform(generator);
			double energy = exponential(generator);
			double amplitude = sqrt(energy);
			chCoef[nuser][i][0] = amplitude * cos(theta); // real part
			chCoef[nuser][i][1] = amplitude * sin(theta); // imaginary part
		}
	}
	
	/*chCoef[0][0][0] = 0.219721;
	chCoef[0][0][1] = -0.492365;
	chCoef[0][1][0] = -0.533266;
	chCoef[0][1][1] = 0.420017;
	*/
}

void MultipleAccessChannel(double stdDev, double ***chCoef, double ***tx, double **rx)
{
	for (int i = 0; i < L; i++)
	{
		rx[i][0] = stdDev * normal(generator) / 2; // real
		rx[i][1] = stdDev * normal(generator) / 2; // imaginary
		//rx[i][0] = 0;
		//rx[i][1] = 0;
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			for (int j = 0; j < NUM_TX; j++)
			{
				rx[i][0] += tx[nuser][i][j] * chCoef[nuser][j][0];
				rx[i][1] += tx[nuser][i][j] * chCoef[nuser][j][1];
			}
		}
	}
	/*
	for (int i = 0; i < L; i++)
		cout << rx[i][0] << " ";
	cout << endl;
	for (int i = 0; i < L; i++)
		cout << rx[i][1] << " ";
	cout << endl;
	system("pause");
	*/
}