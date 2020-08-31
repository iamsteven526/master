#define _USE_MATH_DEFINES
#include <iostream>
#include "parameters.h"
#include <math.h>
using namespace std;

void MLDT(double variance, double **chCoef, double **rx, double *app, double **appLlr,double **estimate)
{
	if (CE_METHOD)
	{
		vector<vector<int>> Quntanary(NUM_LEVEL, vector<int>(NUM_USER));
		for (int level = 0; level < NUM_LEVEL; level++)
		{
			int reg = level;
			for (int nuser = 0; nuser < NUM_USER; nuser++)
			{
				Quntanary[level][NUM_USER - nuser - 1] = reg % 4;
				reg /= 4;
			}
		}
	
		vector<vector<double>> L(NUM_LEVEL, vector<double>(2));

		for (int level = 0; level < NUM_LEVEL; level++)
		{
			for (int nuser = 0; nuser < NUM_USER; nuser++)
			{
				L[level][0] += cos(M_PI / 4 + (M_PI / 2) * Quntanary[level][nuser]) * chCoef[nuser][0];
				L[level][0] -= sin(M_PI / 4 + (M_PI / 2) * Quntanary[level][nuser]) * chCoef[nuser][1];
				L[level][1] += cos(M_PI / 4 + (M_PI / 2) * Quntanary[level][nuser]) * chCoef[nuser][1];
				L[level][1] += sin(M_PI / 4 + (M_PI / 2) * Quntanary[level][nuser]) * chCoef[nuser][0];
			}
		}

		
		for (int i = 0; i < BLOCK_LEN; i++)
		{
			//---------- APPs ----------
			
			for (int level = 0; level < NUM_LEVEL; level++)
			{
				app[level] = exp(-pow((rx[i][0] - L[level][0]), 2.) / (2. * variance));
			}
			//---------- normalization ----------
			double temp = 0;
			for (int j = 0; j < NUM_LEVEL; j++)
			{
				if (app[j] < NUMERIC_LIMIT) app[j] = NUMERIC_LIMIT;
				temp += app[j];
			}
			for (int j = 0; j < NUM_LEVEL; j++)
			{
				app[j] /= temp;
				if (app[j] < NUMERIC_LIMIT) app[j] = NUMERIC_LIMIT;
			}
			//---------- APPs ----------
			for (int level = 0; level < NUM_LEVEL; level++)
			{
				app[level] *= exp(-pow((rx[i][1] - L[level][1]), 2.) / (2. * variance));
			}
			//---------- normalization ----------
			temp = 0;
			for (int j = 0; j < NUM_LEVEL; j++)
			{
				if (app[j] < NUMERIC_LIMIT) app[j] = NUMERIC_LIMIT;
				temp += app[j];
			}
			for (int j = 0; j < NUM_LEVEL; j++)
			{
				app[j] /= temp;
				if (app[j] < NUMERIC_LIMIT) app[j] = NUMERIC_LIMIT;
			}
			//---------- a-posteriori LLRs ----------
			for (int nuser = 0; nuser < NUM_USER; nuser++)
			{
				double upper[2] = { 0 }, lower[2] = { 0 };
				for (int level = 0; level < NUM_LEVEL; level++)
				{
					if (Quntanary[level][nuser] == 0)
					{
						upper[0] += app[level];
						upper[1] += app[level];
					}
					else if (Quntanary[level][nuser] == 1)
					{
						lower[0] += app[level];
						upper[1] += app[level];
					}
					else if (Quntanary[level][nuser] == 2)
					{
						lower[0] += app[level];
						lower[1] += app[level];
					}
					else if (Quntanary[level][nuser] == 3)
					{
						upper[0] += app[level];
						lower[1] += app[level];
					}
					else
					{
						cout << "error!";
						system("pause");
					}
				}

				for (int m = 0; m < 2; m++)
				{
					if (upper[m] <=NUMERIC_LIMIT)  appLlr[nuser][2 * i + m] = -LLR_LIMIT;
					else if (lower[m] <= NUMERIC_LIMIT) appLlr[nuser][2 * i + m] = LLR_LIMIT;
					else appLlr[nuser][2 * i + m] = log(upper[m] / lower[m]);
				}
			}
			/*
			if (NUM_USER == 1)
			{
				if ((app[0] + app[1]) <= NUMERIC_LIMIT) appLlr[0][2 * i] = -LLR_LIMIT;
				else if ((app[2] + app[3]) <= NUMERIC_LIMIT) appLlr[0][2 * i] = LLR_LIMIT;
				else appLlr[0][2 * i] = log((app[0] + app[1]) / (app[2] + app[3]));
				if ((app[0] + app[2]) <= NUMERIC_LIMIT) appLlr[0][2 * i + 1] = -LLR_LIMIT;
				else if ((app[1] + app[3]) <= NUMERIC_LIMIT) appLlr[0][2 * i + 1] = LLR_LIMIT;
				else appLlr[0][2 * i + 1] = log((app[0] + app[2]) / (app[1] + app[3]));
			}
			else if (NUM_USER == 2)
			{
				//---------- user-A ----------
				if ((app[0] + app[1] + app[2] + app[3] + app[4] + app[5] + app[6] + app[7]) <= NUMERIC_LIMIT) appLlr[0][2 * i] = -LLR_LIMIT;
				else if ((app[8] + app[9] + app[10] + app[11] + app[12] + app[13] + app[14] + app[15]) <= NUMERIC_LIMIT) appLlr[0][2 * i] = LLR_LIMIT;
				else appLlr[0][2 * i] = log((app[0] + app[1] + app[2] + app[3] + app[4] + app[5] + app[6] + app[7]) / (app[8] + app[9] + app[10] + app[11] + app[12] + app[13] + app[14] + app[15]));
				if ((app[0] + app[1] + app[2] + app[3] + app[8] + app[9] + app[10] + app[11]) <= NUMERIC_LIMIT) appLlr[0][2 * i + 1] = -LLR_LIMIT;
				else if ((app[4] + app[5] + app[6] + app[7] + app[12] + app[13] + app[14] + app[15]) <= NUMERIC_LIMIT) appLlr[0][2 * i + 1] = LLR_LIMIT;
				else appLlr[0][2 * i + 1] = log((app[0] + app[1] + app[2] + app[3] + app[8] + app[9] + app[10] + app[11]) / (app[4] + app[5] + app[6] + app[7] + app[12] + app[13] + app[14] + app[15]));
				//---------- user-B ----------
				if ((app[0] + app[1] + app[4] + app[5] + app[8] + app[9] + app[12] + app[13]) <= NUMERIC_LIMIT) appLlr[1][2 * i] = -LLR_LIMIT;
				else if ((app[2] + app[3] + app[6] + app[7] + app[10] + app[11] + app[14] + app[15]) <= NUMERIC_LIMIT) appLlr[1][2 * i] = LLR_LIMIT;
				else appLlr[1][2 * i] = log((app[0] + app[1] + app[4] + app[5] + app[8] + app[9] + app[12] + app[13]) / (app[2] + app[3] + app[6] + app[7] + app[10] + app[11] + app[14] + app[15]));
				if ((app[0] + app[2] + app[4] + app[6] + app[8] + app[10] + app[12] + app[14]) <= NUMERIC_LIMIT) appLlr[1][2 * i + 1] = -LLR_LIMIT;
				else if ((app[1] + app[3] + app[5] + app[7] + app[9] + app[11] + app[13] + app[15]) <= NUMERIC_LIMIT) appLlr[1][2 * i + 1] = LLR_LIMIT;
				else appLlr[1][2 * i + 1] = log((app[0] + app[2] + app[4] + app[6] + app[8] + app[10] + app[12] + app[14]) / (app[1] + app[3] + app[5] + app[7] + app[9] + app[11] + app[13] + app[15]));
			}
			*/
		}
	}
	else
	{
		for (int i = 0; i < NUM_USER; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				chCoef[i][j]=estimate[i][j];
			}
		}
		for (int i = 0; i < BLOCK_LEN; i++)
		{
			//---------- APPs ----------
			if (NUM_USER == 1)
			{
				app[0] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1])), 2.) / (2. * variance));
				app[1] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1])), 2.) / (2. * variance));
				app[2] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1])), 2.) / (2. * variance));
				app[3] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1])), 2.) / (2. * variance));
			}
			else if (NUM_USER == 2)
			{
				app[0] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[1] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[2] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[3] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[4] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[5] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[6] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[7] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[8] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[9] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[10] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[11] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[12] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[13] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[14] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[15] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
			}
			else
			{
				printf("\nPARAMETER SETTING IS WRONG\n");
				system("pause");
			}
			//---------- normalization ----------
			double temp = 0;
			for (int j = 0; j < NUM_LEVEL; j++)
			{
				if (app[j] < NUMERIC_LIMIT) app[j] = NUMERIC_LIMIT;
				temp += app[j];
			}
			for (int j = 0; j < NUM_LEVEL; j++)
			{
				app[j] /= temp;
				if (app[j] < NUMERIC_LIMIT) app[j] = NUMERIC_LIMIT;
			}
			//---------- APPs ----------
			if (NUM_USER == 1)
			{
				app[0] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1])), 2.) / (2. * variance));
				app[1] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1])), 2.) / (2. * variance));
				app[2] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1])), 2.) / (2. * variance));
				app[3] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1])), 2.) / (2. * variance));
			}
			else if (NUM_USER == 2)
			{
				app[0] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[1] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[2] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[3] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[4] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[5] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[6] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[7] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[8] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[9] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[10] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[11] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[12] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[13] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[14] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[15] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
			}
			//---------- normalization ----------
			temp = 0;
			for (int j = 0; j < NUM_LEVEL; j++)
			{
				if (app[j] < NUMERIC_LIMIT) app[j] = NUMERIC_LIMIT;
				temp += app[j];
			}
			for (int j = 0; j < NUM_LEVEL; j++)
			{
				app[j] /= temp;
				if (app[j] < NUMERIC_LIMIT) app[j] = NUMERIC_LIMIT;
			}
			//---------- a-posteriori LLRs ----------
			if (NUM_USER == 1)
			{
				if ((app[0] + app[1]) <= NUMERIC_LIMIT) appLlr[0][2 * i] = -LLR_LIMIT;
				else if ((app[2] + app[3]) <= NUMERIC_LIMIT) appLlr[0][2 * i] = LLR_LIMIT;
				else appLlr[0][2 * i] = log((app[0] + app[1]) / (app[2] + app[3]));
				if ((app[0] + app[2]) <= NUMERIC_LIMIT) appLlr[0][2 * i + 1] = -LLR_LIMIT;
				else if ((app[1] + app[3]) <= NUMERIC_LIMIT) appLlr[0][2 * i + 1] = LLR_LIMIT;
				else appLlr[0][2 * i + 1] = log((app[0] + app[2]) / (app[1] + app[3]));
			}
			else if (NUM_USER == 2)
			{
				//---------- user-A ----------
				if ((app[0] + app[1] + app[2] + app[3] + app[4] + app[5] + app[6] + app[7]) <= NUMERIC_LIMIT) appLlr[0][2 * i] = -LLR_LIMIT;
				else if ((app[8] + app[9] + app[10] + app[11] + app[12] + app[13] + app[14] + app[15]) <= NUMERIC_LIMIT) appLlr[0][2 * i] = LLR_LIMIT;
				else appLlr[0][2 * i] = log((app[0] + app[1] + app[2] + app[3] + app[4] + app[5] + app[6] + app[7]) / (app[8] + app[9] + app[10] + app[11] + app[12] + app[13] + app[14] + app[15]));
				if ((app[0] + app[1] + app[2] + app[3] + app[8] + app[9] + app[10] + app[11]) <= NUMERIC_LIMIT) appLlr[0][2 * i + 1] = -LLR_LIMIT;
				else if ((app[4] + app[5] + app[6] + app[7] + app[12] + app[13] + app[14] + app[15]) <= NUMERIC_LIMIT) appLlr[0][2 * i + 1] = LLR_LIMIT;
				else appLlr[0][2 * i + 1] = log((app[0] + app[1] + app[2] + app[3] + app[8] + app[9] + app[10] + app[11]) / (app[4] + app[5] + app[6] + app[7] + app[12] + app[13] + app[14] + app[15]));
				//---------- user-B ----------
				if ((app[0] + app[1] + app[4] + app[5] + app[8] + app[9] + app[12] + app[13]) <= NUMERIC_LIMIT) appLlr[1][2 * i] = -LLR_LIMIT;
				else if ((app[2] + app[3] + app[6] + app[7] + app[10] + app[11] + app[14] + app[15]) <= NUMERIC_LIMIT) appLlr[1][2 * i] = LLR_LIMIT;
				else appLlr[1][2 * i] = log((app[0] + app[1] + app[4] + app[5] + app[8] + app[9] + app[12] + app[13]) / (app[2] + app[3] + app[6] + app[7] + app[10] + app[11] + app[14] + app[15]));
				if ((app[0] + app[2] + app[4] + app[6] + app[8] + app[10] + app[12] + app[14]) <= NUMERIC_LIMIT) appLlr[1][2 * i + 1] = -LLR_LIMIT;
				else if ((app[1] + app[3] + app[5] + app[7] + app[9] + app[11] + app[13] + app[15]) <= NUMERIC_LIMIT) appLlr[1][2 * i + 1] = LLR_LIMIT;
				else appLlr[1][2 * i + 1] = log((app[0] + app[2] + app[4] + app[6] + app[8] + app[10] + app[12] + app[14]) / (app[1] + app[3] + app[5] + app[7] + app[9] + app[11] + app[13] + app[15]));
			}
		}
	}
}