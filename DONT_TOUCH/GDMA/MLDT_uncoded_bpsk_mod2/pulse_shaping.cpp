#define _USE_MATH_DEFINES
#include <iostream>
#include "parameters.h"
#include <string>
#include <fstream>
#include <cmath>
#include <cstring>
using namespace std;
#pragma warning(disable:4996)
void UpSampling(double *tx, double *preTx, double *polit, double *prePilot)
{
	memset(preTx, 0, sizeof(double)*BLOCK_LEN*UP_RATE*SPREAD_LEN);
	for (int i = 0; i < BLOCK_LEN * SPREAD_LEN; i++)
	{
		preTx[i*UP_RATE] = tx[i];
		prePilot[i * UP_RATE] = polit[i];
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
	cout << drift << " ";
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
	//print_shape(srrc);
}

void PulseShaping(double *preTx, double *tx, double *txFilter,double *rxFilter,double *prePilot,double *pilot)
{
	int effLen = BLOCK_LEN*UP_RATE*SPREAD_LEN;
	int trunLen = TRUNCATION*UP_RATE;
	memset(tx, 0, sizeof(double) * effLen);
	for (int i = 0; i < effLen; i++) // transmitting filter
	{
		for (int j = -trunLen; j < trunLen + UP_RATE; j++) // convolution
		{
			tx[i] += txFilter[j + trunLen] * preTx[(effLen + (i - j) % effLen) % effLen];
		}
	}
	
	memset(pilot, 0, sizeof(double) * effLen);
	for (int i = 0; i < effLen; i++) // transmitting filter
	{
		for (int j = -trunLen; j < trunLen + UP_RATE; j++) // convolution
		{
			pilot[i] += txFilter[j + trunLen] * prePilot[(effLen + (i - j) % effLen) % effLen];
		}
	}
//	print_shape(tx,"syn");
	/*memset(tx, 0, sizeof(double)*effLen);
	for (int i = 0; i < effLen; i++) // transmitting filter
	{
		for (int j = -trunLen; j < trunLen + UP_RATE; j++) // convolution
		{
			tx[i] += txFilter[j + trunLen] * preTx[(effLen + (i - j) % effLen) % effLen];
		}
	}*/
//	print_shape(tx,"drift");

}

void ReceivingFilter(double **rx, double **postRx, double *rxFilter, double *txFilter, double** pilot,double **chCoef)
{

	int effLen = BLOCK_LEN*UP_RATE*SPREAD_LEN;
	int trunLen = TRUNCATION*UP_RATE;
	double estimate_drift[NUM_USER] = { 0 }, m[NUM_USER] = { 0 };

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
		/*for (int j = 0; j < BLOCK_LEN; j++) // down-sampling
		{
			rx[j][i] = postRx[j*UP_RATE][i];
		}*/
	}

	//-------------m_k(0) calculation

	double max = 0;
	double temp = 0;
	double mean[NUM_USER][UP_RATE][2] = { 0 };
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		for (int i = 0; i < UP_RATE; i++)
		{
			for (int j = 0; j < SAMPLE_NUM; j++)
			{
				temp += pow(postRx[j * UP_RATE + i][0] * pilot[nuser][j], 2) + pow(postRx[j * UP_RATE + i][1] * pilot[nuser][j], 2);
			}
			
			//cout << temp << " ";
			if (temp > max)
			{
				m[nuser] = i;
				max = temp;
			}
			temp = 0;
		}
		estimate_drift[nuser] = m[nuser] / UP_RATE;
		//cout << endl;
		//cout << m[nuser];
		cout << estimate_drift[nuser] << " ";
	}

	//--Delayed Locked Loop
	double a[NUM_USER][UP_RATE] = { 0 }, b[NUM_USER][UP_RATE] = { 0 };
	
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		for (int i = 0; i < UP_RATE; i++)
		{
			a[nuser][i] = chCoef[nuser][0] * (postRx[i][0] - mean[nuser][i][0]) - chCoef[nuser][1] * (postRx[i][1] - mean[nuser][i][1]);
			b[nuser][i] = chCoef[nuser][0] * (postRx[i][1] - mean[nuser][i][1]) + chCoef[nuser][1] * (postRx[i][0] - mean[nuser][i][0]);
		}
	}




	//system("pause");
	/*ofstream syn_real("syn_real.txt");

	for (int i = 0; i < 20 * UP_RATE; i++)
		syn_real << postRx[i][0] << " ";*/
	



/*	ofstream drift_real("drift_real.txt");
	for (int i = 0; i < 20* UP_RATE; i++)
		drift_real << postRx[i][0] << " ";

	ofstream drift_img("drift_img.txt");
	for (int i = 0; i < 20 * UP_RATE; i++)
		drift_img << postRx[i][1] << " ";*/


/*	for (int i = 0; i < 2; i++) // real and imaginary
	{
		for (int j = 0; j < effLen; j++) // receiving filter
		{
			postRx[j][i] = 0;
			for (int k = -trunLen; k < trunLen + UP_RATE; k++)
			{
				postRx[j][i] += txFilter[k + trunLen] * rx[(effLen + (j - k) % effLen) % effLen][i];
			}
		}
		for (int j = 0; j < BLOCK_LEN; j++) // down-sampling
		{
			rx[j][i] = postRx[j * UP_RATE][i];
		}
	}*/

/*	ofstream syn_real("syn_real.txt");
	for (int i = 0; i < 20 * UP_RATE; i++)
		syn_real << postRx[i][0] << " ";

	ofstream syn_img("syn_img.txt");
	for (int i = 0; i < 20 * UP_RATE; i++)
		syn_img << postRx[i][1] << " ";

	system("pause");*/
}

void print_shape(double* txfilter, string name)
{
	string ext = ".txt";
	name += ext;
	ofstream outfile(name.c_str());

	for (int i=0;i<20* UP_RATE;i++)
	{
		outfile<< txfilter[i]<<" ";
	}
}