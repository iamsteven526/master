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
	if (Racian == 0)
	{
		if (FADING_TYPE == 0 || FADING_TYPE == 2)
		{
			//---- block fading
			for (int i = 0; i < int(CODE_LEN) / int(FADING_SIZE); i++)
			{
				for (int nuser = 0; nuser < NUM_USER; nuser++)
				{
					double theta = 2 * M_PI * uniform(generator);		// phase shift
					chCoef[nuser][3][i] = exponential(generator);			// energy
					//chCoef[nuser][3][i] = 1;
					chCoef[nuser][2][i] = sqrt(chCoef[nuser][3][i]);			// amplitide
					chCoef[nuser][1][i] = chCoef[nuser][2][i] * sin(theta);	// imaginary pary
					chCoef[nuser][0][i] = chCoef[nuser][2][i] * cos(theta);	// real part
				}
			}
			/*chCoef[0][3][0] = 1;
			chCoef[0][2][0] = sqrt(chCoef[0][3][0]);			// amplitide
			chCoef[0][1][0] = chCoef[0][2][0] * sin(0);	// imaginary pary
			chCoef[0][0][0] = chCoef[0][2][0] * cos(0);	// real part
			
			chCoef[1][3][0] = 1;
			chCoef[1][2][0] = sqrt(chCoef[1][3][0]);			// amplitide
			chCoef[1][1][0] = chCoef[1][2][0] * sin(M_PI / 2);	// imaginary pary
			chCoef[1][0][0] = chCoef[1][2][0] * cos(M_PI / 2);	// real part*/

		//	system("pause");
		}
		else
		{
			//---- Doppler fading
			double period = 1. / (double)(CHIP_RATE * pow(10, 6)); // sampling duration
			double fd = (UE_SPEED * 0.277778) * (CARRIER_FREQ * pow(10, 9)) / (3. * pow(10, 8)); // maximun Doppler frequency shift
			const int M = 16;
			double nFactor = 0;

			for (int i = 0; i < NUM_USER; i++)
			{

				double theta = uniform(generator) * 2. * M_PI - M_PI;
				for (int n = 1; n <= M; n++)
				{
					double phi_1 = uniform(generator) * 2. * M_PI - M_PI;
					double phi_2 = uniform(generator) * 2. * M_PI - M_PI;
					double alpha = (2. * M_PI * n - M_PI + theta) / (4. * M);
					for (int k = 0; k < CODE_LEN; k++)
					{
						double t = k * period;
						chCoef[i][0][k] += cos(2. * M_PI * fd * t * cos(alpha) + phi_1);
						chCoef[i][1][k] += cos(2. * M_PI * fd * t * sin(alpha) + phi_1);
					}
				}

				// Derive different channel tap
				double phi = uniform(generator) * 2. * M_PI - M_PI;
				double thi = uniform(generator) * 2. * M_PI - M_PI;
				for (int k = 0; k <CODE_LEN; k++)
				{
					double t = k * period;
					chCoef[i][0][k] = pow(K + 1, -0.5) * chCoef[i][0][k] * sqrt(2. / M) + pow(K, 0.5) / pow(K + 1, 0.5) * cos(2. * M_PI * fd * t * cos(thi) + phi);
					chCoef[i][1][k] = pow(K + 1, -0.5) * chCoef[i][1][k] * sqrt(2. / M) + pow(K, 0.5) / pow(K + 1, 0.5) * sin(2. * M_PI * fd * t * cos(thi) + phi);
				}
			}
		}
	}
	else 
	{
		if (NUM_USER != 2)
		{
			double phi = 2 * M_PI * uniform(generator);
			double theta = 2 * M_PI * uniform(generator);																	// phase shift

			chCoef[0][3][0] = exponential(generator);																		// energy
			chCoef[0][2][0] = sqrt(chCoef[0][3][0]);																		// amplitide
			chCoef[0][1][0] = pow(K, 0.5) / pow(K + 1, 0.5) * sin(phi) + chCoef[0][2][0] * sin(theta) / pow(K + 1, 0.5);	// imaginary pary
			chCoef[0][0][0] = pow(K, 0.5) / pow(K + 1, 0.5) * cos(phi) + chCoef[0][2][0] * cos(theta) / pow(K + 1, 0.5);	// real part
			
		}
		else
		{
			double phi_1, phi_2;

			if (phase != 11)
			{
				phi_1 = 0;
				phi_2 = M_PI * phase / 180;
			}
			else
			{
				phi_1 = 2 * M_PI * uniform(generator);
				phi_2 = 2 * M_PI * uniform(generator);
			}

			for (int i = 0; i < CODE_LEN / FADING_SIZE + 1; i++)
			{
				double theta = 2 * M_PI * uniform(generator);																	// phase shift

				chCoef[0][3][i] = exponential(generator);																		// energy
				chCoef[0][2][i] = sqrt(chCoef[0][3][i]);																		// amplitide
				chCoef[0][1][i] = pow(K, 0.5) / pow(K + 1, 0.5) * sin(phi_1) + chCoef[0][2][i] * sin(theta) / pow(K + 1, 0.5);	// imaginary pary
				chCoef[0][0][i] = pow(K, 0.5) / pow(K + 1, 0.5) * cos(phi_1) + chCoef[0][2][i] * cos(theta) / pow(K + 1, 0.5);	// real part

				theta = 2 * M_PI * uniform(generator);																			// phase shift
				chCoef[1][3][i] = exponential(generator);																		// energy
				chCoef[1][2][i] = sqrt(chCoef[0][3][i]);																		// amplitide
				chCoef[1][1][i] = pow(K, 0.5) / pow(K + 1, 0.5) * sin(phi_2) + chCoef[0][2][i] * sin(theta) / pow(K + 1, 0.5);	// imaginary pary
				chCoef[1][0][i] = pow(K, 0.5) / pow(K + 1, 0.5) * cos(phi_2) + chCoef[0][2][i] * cos(theta) / pow(K + 1, 0.5);	// real part	

			}
		}
		
		//block fading
		/*for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			for (int i = 1; i < CODE_LEN; i++)
			{
				chCoef[nuser][1][i] = chCoef[nuser][1][0];
				chCoef[nuser][0][i] = chCoef[nuser][0][0];
			}
		}*/
	}
}

void MultipleAccessChannel(double stdDev, double ***chCoef, double **tx, double **rx)
{
	int effLen = SYNCHRONOUS ? CODE_LEN : CODE_LEN*UP_RATE;

	for (int i = 0; i < effLen; i++)
	{
		rx[i][0] = stdDev * normal(generator); // real
		rx[i][1] = stdDev * normal(generator); // imaginary
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			rx[i][0] += tx[nuser][i] * chCoef[nuser][0][i / FADING_SIZE];
			rx[i][1] += tx[nuser][i] * chCoef[nuser][1][i / FADING_SIZE];
		}
	}

}
