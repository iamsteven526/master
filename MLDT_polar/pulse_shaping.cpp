#define _USE_MATH_DEFINES
#include <iostream>
#include "parameters.h"
using namespace std;

void UpSampling(double *tx, double *preTx)
{
	memset(preTx, 0, sizeof(double)*CODE_LEN*UP_RATE);
	for (int i = 0; i < CODE_LEN; i++)
	{
		preTx[i*UP_RATE] = tx[i];
	}
}

double SquareRootRaisedCosine(double m)
{
	if (m == 0)
	{
		return (1. + ROLL_OFF*((4. / M_PI) - 1.));
	}
	else if (fabs(m - (1. / (4.*ROLL_OFF))) < 1e-8 || fabs(m + (1. / (4.*ROLL_OFF))) < 1e-8)
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
			double m = (double)j / (double)UP_RATE + n;
			srrc[i * UP_RATE + j] = (fabs(m) < TRUNCATION) ? SquareRootRaisedCosine(m) : 0;
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
	}
}

void PulseShaping(double *preTx, double *tx, double *txFilter)
{
	int effLen = CODE_LEN*UP_RATE;
	int trunLen = TRUNCATION*UP_RATE;
	memset(tx, 0, sizeof(double)*effLen);
	for (int i = 0; i < effLen; i++) // transmitting filter
	{
		for (int j = -trunLen; j < trunLen + UP_RATE; j++)
		{
			tx[i] += txFilter[j + trunLen] * preTx[(effLen + (i - j) % effLen) % effLen];
		}
	}
}

void ReceivingFilter(double **rx, double **postRx, double *rxFilter)
{
	int effLen = CODE_LEN*UP_RATE;
	int trunLen = TRUNCATION*UP_RATE;
	for (int i = 0; i < 2; i++) // real and imaginary
	{
		for (int j = 0; j < effLen; j++) // receiving filter
		{
			postRx[j][i] = 0;
			for (int k = -trunLen; k < trunLen + UP_RATE; k++)
			{
				postRx[j][i] += rxFilter[k + trunLen] * rx[(effLen + (j - k) % effLen) % effLen][i];
			}
		}
		for (int j = 0; j < CODE_LEN; j++) // down-sampling
		{
			rx[j][i] = postRx[j*UP_RATE][i];
		}
	}
}