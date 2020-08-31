#define _USE_MATH_DEFINES
#include <iostream>
#include <random>
#include "parameters.h"
#include <cmath>
#include <cstring>
using namespace std;

namespace
{
	random_device seed;
	mt19937 generator(seed());
	uniform_real_distribution<double> uniform(-0.5*DRIFT_RANGE, 0.5*DRIFT_RANGE);
}

void UpSampling(double ****tx, double ****preTx)
{
	if (!OVER_SAMPLE)	 //----------------����over sampling�A�u���Ҳv pulse shaping + CP or CS 
	{
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			for (int i = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
			{
				if (CP_TYPE)			 //---------------- pulse shaping + CP
				{
					for (int j = 0; j < (CP_LEN + FFT_POINT) * UP_RATE; j++)
					{
						preTx[nuser][i][j][0] = preTx[nuser][i][j][1] = 0;
					}
					for (int j = 0; j < CP_LEN + FFT_POINT; j++)
					{
						preTx[nuser][i][j * UP_RATE][0] = tx[nuser][i][j][0];
						preTx[nuser][i][j * UP_RATE][1] = tx[nuser][i][j][1];
					}
				}
				else					//---------------- pulse shaping + CP + CS
				{
					for (int j = 0; j < (CP_LEN + CS_LEN + FFT_POINT) * UP_RATE; j++)
					{
						preTx[nuser][i][j][0] = preTx[nuser][i][j][1] = 0;
					}
					for (int j = 0; j < CP_LEN + CS_LEN + FFT_POINT; j++)
					{
						preTx[nuser][i][j * UP_RATE][0] = tx[nuser][i][j][0];
						preTx[nuser][i][j * UP_RATE][1] = tx[nuser][i][j][1];
					}
				}
			}
		}
	}
	else
	{
		if (PULSE_SHAPE)				
		{
			for (int nuser = 0; nuser < NUM_USER; nuser++)
			{
				for (int i = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
				{
					if (CP_TYPE)
					{
						for (int j = 0; j < (CP_LEN + FFT_POINT) * OVER_SAMPLE_RATE * UP_RATE; j++)
						{
							preTx[nuser][i][j][0] = preTx[nuser][i][j][1] = 0;
						}
						for (int j = 0; j < (CP_LEN + FFT_POINT) * OVER_SAMPLE_RATE; j++)
						{
							preTx[nuser][i][j * UP_RATE][0] = tx[nuser][i][j][0];
							preTx[nuser][i][j * UP_RATE][1] = tx[nuser][i][j][1];
						}
					}
					else		//----------------over_sampling + pulse shaping + CP + CS
					{
						for (int j = 0; j < (CP_LEN + CS_LEN + FFT_POINT) * OVER_SAMPLE_RATE * UP_RATE; j++)
						{
							preTx[nuser][i][j][0] = preTx[nuser][i][j][1] = 0;
						}
						for (int j = 0; j < (CP_LEN + CS_LEN + FFT_POINT) * OVER_SAMPLE_RATE; j++)
						{
							preTx[nuser][i][j * UP_RATE][0] = tx[nuser][i][j][0];
							preTx[nuser][i][j * UP_RATE][1] = tx[nuser][i][j][1];
						}

					}
				}
			}
		}
		else					//----------------over_sampling + CP + CS
		{
			for (int nuser = 0; nuser < NUM_USER; nuser++)
			{
				for (int i = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
				{
					if (CP_TYPE)
					{
						for (int j = 0; j < (CP_LEN + FFT_POINT) * OVER_SAMPLE_RATE; j++)
						{
							preTx[nuser][i][j][0] = tx[nuser][i][j][0];
							preTx[nuser][i][j][1] = tx[nuser][i][j][1];
						}
					}
					else
					{
						for (int j = 0; j < (CP_LEN + CS_LEN + FFT_POINT) * OVER_SAMPLE_RATE; j++)
						{
							preTx[nuser][i][j][0] = tx[nuser][i][j][0];
							preTx[nuser][i][j][1] = tx[nuser][i][j][1];
						}
					}
				}
			}
		}
	}
}

double SquareRootRaisedCosine(double m)
{
	if (m == 0)
	{
		return (1. + ROLL_OFF*(4. / M_PI - 1.));
	}
	else if (fabs(m - 1. / (4.*ROLL_OFF)) < 1e-8 || fabs(m + 1. / (4.*ROLL_OFF)) < 1e-8)
	{
		return ((ROLL_OFF / sqrt(2.))*((1. + 2. / M_PI)*sin(M_PI / 4. / ROLL_OFF) + (1. - 2. / M_PI)*cos(M_PI / 4. / ROLL_OFF)));
	}
	else
	{
		return ((sin((1. - ROLL_OFF)*M_PI*m) + 4.*ROLL_OFF*m*cos((1. + ROLL_OFF)*M_PI*m)) / (M_PI*m*(1. - pow(4.*ROLL_OFF*m, 2.))));
	}
}

void SRRCGeneration(double *srrc, double drift)
{
	for (int i = 0; i <= 2 * TRUNCATION; i++)
	{
		double n = i - TRUNCATION - drift;
		for (int j = 0; j < UP_RATE; j++)
		{
			double m = (double)j / UP_RATE + n;
			srrc[i * UP_RATE + j] = fabs(m) < TRUNCATION ? SquareRootRaisedCosine(m) : 0;
		}
	}
	//---------- normalization ---------
	static bool flag = true;
	static double nFactor;
	if (flag)
	{
		nFactor = 0;
		for (int i = 0; i < (2 * TRUNCATION + 1)*UP_RATE; i++)
		{
			nFactor += pow(srrc[i], 2);
		}
		nFactor = sqrt(nFactor);
		flag = false;
	}
	for (int i = 0; i < (2 * TRUNCATION + 1)*UP_RATE; i++)
	{
		srrc[i] /= nFactor;
		//if (OVER_SAMPLE)
		//	srrc[i] *= sqrt(1/(double)OVER_SAMPLE_RATE);
	}
}

void PulseShaping(double ****preTx, double ****tx, double *txFilter, double* drift)
{
	int effLen = OVER_SAMPLE ? CP_TYPE ? (CP_LEN + FFT_POINT) * OVER_SAMPLE_RATE *  UP_RATE : (CP_LEN + CS_LEN + FFT_POINT) * OVER_SAMPLE_RATE * UP_RATE : CP_TYPE ? (CP_LEN + FFT_POINT) * UP_RATE : (CP_LEN + CS_LEN + FFT_POINT) * UP_RATE;
	int trunLen = TRUNCATION * UP_RATE;
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		//drift[nuser] = uniform(generator);
		SRRCGeneration(txFilter, 0);
		for (int i = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
		{
			for (int j = 0; j < effLen; j++) // transmitting filter
			{
				tx[nuser][i][j][0] = tx[nuser][i][j][1] = 0;
				for (int k = -trunLen; k < trunLen + UP_RATE; k++)
				{
					tx[nuser][i][j][0] += txFilter[k + trunLen] * preTx[nuser][i][(effLen + (j - k) % effLen) % effLen][0]; // real
					tx[nuser][i][j][1] += txFilter[k + trunLen] * preTx[nuser][i][(effLen + (j - k) % effLen) % effLen][1]; // imaginary
				}
			}
		}
	}

}

void ReceivingFilter(double ***rx, double **tempRx, double *rxFilter)
{
	int effLen = OVER_SAMPLE ? CP_TYPE ? (CP_LEN + FFT_POINT) * OVER_SAMPLE_RATE * UP_RATE : (CP_LEN + CS_LEN + FFT_POINT) * OVER_SAMPLE_RATE * UP_RATE : CP_TYPE ? (CP_LEN + FFT_POINT) * UP_RATE : (CP_LEN + CS_LEN + FFT_POINT) * UP_RATE;
	
	int trunLen = TRUNCATION * UP_RATE;

	for (int i = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
	{
		for (int j = 0; j < effLen; j++) // receiving filter
		{
			tempRx[j][0] = tempRx[j][1] = 0;
			for (int k = -trunLen; k < trunLen + UP_RATE; k++)
			{
				tempRx[j][0] += rxFilter[k + trunLen] * rx[i][(effLen + (j - k) % effLen) % effLen][0]; // real
				tempRx[j][1] += rxFilter[k + trunLen] * rx[i][(effLen + (j - k) % effLen) % effLen][1]; // imaginary
			}
		}
			
		if (!OVER_SAMPLE)
		{
			if (CP_TYPE)
			{
				for (int j = 0; j < FFT_POINT + CP_LEN; j++) // down-sampling
				{
					rx[i][j][0] = tempRx[j * UP_RATE][0];
					rx[i][j][1] = tempRx[j * UP_RATE][1];
				}
			}
			else
			{
				for (int j = 0; j < FFT_POINT + CP_LEN + CS_LEN; j++) // down-sampling
				{
					rx[i][j][0] = tempRx[j * UP_RATE][0];
					rx[i][j][1] = tempRx[j * UP_RATE][1];
				}
			}
		}
		else
		{
			if (CP_TYPE)
			{
				for (int j = 0; j < (FFT_POINT + CP_LEN) * OVER_SAMPLE_RATE; j++) // down-sampling
				{
					rx[i][j][0] = tempRx[j * UP_RATE][0];
					rx[i][j][1] = tempRx[j * UP_RATE][1];
				}
			}
			else
			{
				for (int j = 0; j < (FFT_POINT + CP_LEN + CS_LEN) * OVER_SAMPLE_RATE; j++) // down-sampling
				{
					rx[i][j][0] = tempRx[j * UP_RATE][0];
					rx[i][j][1] = tempRx[j * UP_RATE][1];
				}
			}
		}
	}
	

}