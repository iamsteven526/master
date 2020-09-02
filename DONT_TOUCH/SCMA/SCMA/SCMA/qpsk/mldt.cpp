#define _USE_MATH_DEFINES
#include <iostream>
#include "parameters.h"
using namespace std;

void MLDT(double variance, double **chCoef, double **rx, double *app, double **appLlr,double **estimate)
{
	if (CE_METHOD)
	{
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