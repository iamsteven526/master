#define _USE_MATH_DEFINES
#include <iostream>
#include "parameters.h"
using namespace std;

void MLDT(double variance, double **chCoef, double **rx, double *app, double **appLlr, double **estimate)
{
	int count = 0;
	for (int i = 0; i < BLOCK_LEN; i++)
	{
		while (true)
		{
			bool flag = true;
			//---------- APPs ----------
			if (NUM_USER == 1)
			{
				if (CE_SCHEME == 1)
				{
					app[0] = exp(-pow(rx[i][0] - chCoef[0][0], 2) / (2. * (variance + count*MODIFIER)));
					app[1] = exp(-pow(rx[i][0] + chCoef[0][0], 2) / (2. * (variance + count*MODIFIER)));
				}
				else
				{
					app[0] = exp(-pow(rx[i][0] - estimate[0][0], 2) / (2. * (variance + count*MODIFIER)));
					app[1] = exp(-pow(rx[i][0] + estimate[0][0], 2) / (2. * (variance + count*MODIFIER)));
				}
			}
			else if (NUM_USER == 2)
			{
				if (CE_SCHEME == 1)
				{
					app[0] = exp(-pow(rx[i][0] - chCoef[0][0] - chCoef[1][0], 2) / (2. * (variance + count*MODIFIER)));
					app[1] = exp(-pow(rx[i][0] - chCoef[0][0] + chCoef[1][0], 2) / (2. * (variance + count*MODIFIER)));
					app[2] = exp(-pow(rx[i][0] + chCoef[0][0] - chCoef[1][0], 2) / (2. * (variance + count*MODIFIER)));
					app[3] = exp(-pow(rx[i][0] + chCoef[0][0] + chCoef[1][0], 2) / (2. * (variance + count*MODIFIER)));
				}
				else
				{
					app[0] = exp(-pow(rx[i][0] - estimate[0][0] - estimate[1][0], 2) / (2. * (variance + count*MODIFIER)));
					app[1] = exp(-pow(rx[i][0] - estimate[0][0] + estimate[1][0], 2) / (2. * (variance + count*MODIFIER)));
					app[2] = exp(-pow(rx[i][0] + estimate[0][0] - estimate[1][0], 2) / (2. * (variance + count*MODIFIER)));
					app[3] = exp(-pow(rx[i][0] + estimate[0][0] + estimate[1][0], 2) / (2. * (variance + count*MODIFIER)));
				}
			}
			else if (NUM_USER == 3)
			{
				if (CE_SCHEME == 1)
				{
					app[0] = exp(-pow(rx[i][0] - chCoef[0][0] - chCoef[1][0] - chCoef[2][0], 2) / (2. * (variance + count*MODIFIER)));
					app[1] = exp(-pow(rx[i][0] - chCoef[0][0] - chCoef[1][0] + chCoef[2][0], 2) / (2. * (variance + count*MODIFIER)));
					app[2] = exp(-pow(rx[i][0] - chCoef[0][0] + chCoef[1][0] - chCoef[2][0], 2) / (2. * (variance + count*MODIFIER)));
					app[3] = exp(-pow(rx[i][0] - chCoef[0][0] + chCoef[1][0] + chCoef[2][0], 2) / (2. * (variance + count*MODIFIER)));
					app[4] = exp(-pow(rx[i][0] + chCoef[0][0] - chCoef[1][0] - chCoef[2][0], 2) / (2. * (variance + count*MODIFIER)));
					app[5] = exp(-pow(rx[i][0] + chCoef[0][0] - chCoef[1][0] + chCoef[2][0], 2) / (2. * (variance + count*MODIFIER)));
					app[6] = exp(-pow(rx[i][0] + chCoef[0][0] + chCoef[1][0] - chCoef[2][0], 2) / (2. * (variance + count*MODIFIER)));
					app[7] = exp(-pow(rx[i][0] + chCoef[0][0] + chCoef[1][0] + chCoef[2][0], 2) / (2. * (variance + count*MODIFIER)));
				}
				else
				{
					app[0] = exp(-pow(rx[i][0] - estimate[0][0] - estimate[1][0] - estimate[2][0], 2) / (2. * (variance + count*MODIFIER)));
					app[1] = exp(-pow(rx[i][0] - estimate[0][0] - estimate[1][0] + estimate[2][0], 2) / (2. * (variance + count*MODIFIER)));
					app[2] = exp(-pow(rx[i][0] - estimate[0][0] + estimate[1][0] - estimate[2][0], 2) / (2. * (variance + count*MODIFIER)));
					app[3] = exp(-pow(rx[i][0] - estimate[0][0] + estimate[1][0] + estimate[2][0], 2) / (2. * (variance + count*MODIFIER)));
					app[4] = exp(-pow(rx[i][0] + estimate[0][0] - estimate[1][0] - estimate[2][0], 2) / (2. * (variance + count*MODIFIER)));
					app[5] = exp(-pow(rx[i][0] + estimate[0][0] - estimate[1][0] + estimate[2][0], 2) / (2. * (variance + count*MODIFIER)));
					app[6] = exp(-pow(rx[i][0] + estimate[0][0] + estimate[1][0] - estimate[2][0], 2) / (2. * (variance + count*MODIFIER)));
					app[7] = exp(-pow(rx[i][0] + estimate[0][0] + estimate[1][0] + estimate[2][0], 2) / (2. * (variance + count*MODIFIER)));
				}
			}
			else
			{
				printf("\nPARAMETER SETTING IS WRONG\n");
				system("pause");
			}
			if (MODIFIED_LLR)
			{
				for (int j = 0; j < NUM_LEVEL; j++)
				{
					if (app[j] == 0)
					{
						flag = false;
						break;
					}
				}
				if (!flag)
				{
					count++;
					continue;
				}
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
				if (CE_SCHEME == 1)
				{
					app[0] *= exp(-pow(rx[i][1] - chCoef[0][1], 2) / (2. * (variance + count*MODIFIER)));
					app[1] *= exp(-pow(rx[i][1] + chCoef[0][1], 2) / (2. * (variance + count*MODIFIER)));
				}
				else
				{
					app[0] *= exp(-pow(rx[i][1] - estimate[0][1], 2) / (2. * (variance + count*MODIFIER)));
					app[1] *= exp(-pow(rx[i][1] + estimate[0][1], 2) / (2. * (variance + count*MODIFIER)));
				}
			}
			else if (NUM_USER == 2)
			{
				if (CE_SCHEME == 1)
				{
					app[0] *= exp(-pow(rx[i][1] - chCoef[0][1] - chCoef[1][1], 2) / (2. * (variance + count*MODIFIER)));
					app[1] *= exp(-pow(rx[i][1] - chCoef[0][1] + chCoef[1][1], 2) / (2. * (variance + count*MODIFIER)));
					app[2] *= exp(-pow(rx[i][1] + chCoef[0][1] - chCoef[1][1], 2) / (2. * (variance + count*MODIFIER)));
					app[3] *= exp(-pow(rx[i][1] + chCoef[0][1] + chCoef[1][1], 2) / (2. * (variance + count*MODIFIER)));
				}
				else
				{
					app[0] *= exp(-pow(rx[i][1] - estimate[0][1] - estimate[1][1], 2) / (2. * (variance + count*MODIFIER)));
					app[1] *= exp(-pow(rx[i][1] - estimate[0][1] + estimate[1][1], 2) / (2. * (variance + count*MODIFIER)));
					app[2] *= exp(-pow(rx[i][1] + estimate[0][1] - estimate[1][1], 2) / (2. * (variance + count*MODIFIER)));
					app[3] *= exp(-pow(rx[i][1] + estimate[0][1] + estimate[1][1], 2) / (2. * (variance + count*MODIFIER)));
				}
			}
			else if (NUM_USER == 3)
			{
				if (CE_SCHEME == 1)
				{
					app[0] *= exp(-pow(rx[i][1] - chCoef[0][1] - chCoef[1][1] - chCoef[2][1], 2) / (2. * (variance + count*MODIFIER)));
					app[1] *= exp(-pow(rx[i][1] - chCoef[0][1] - chCoef[1][1] + chCoef[2][1], 2) / (2. * (variance + count*MODIFIER)));
					app[2] *= exp(-pow(rx[i][1] - chCoef[0][1] + chCoef[1][1] - chCoef[2][1], 2) / (2. * (variance + count*MODIFIER)));
					app[3] *= exp(-pow(rx[i][1] - chCoef[0][1] + chCoef[1][1] + chCoef[2][1], 2) / (2. * (variance + count*MODIFIER)));
					app[4] *= exp(-pow(rx[i][1] + chCoef[0][1] - chCoef[1][1] - chCoef[2][1], 2) / (2. * (variance + count*MODIFIER)));
					app[5] *= exp(-pow(rx[i][1] + chCoef[0][1] - chCoef[1][1] + chCoef[2][1], 2) / (2. * (variance + count*MODIFIER)));
					app[6] *= exp(-pow(rx[i][1] + chCoef[0][1] + chCoef[1][1] - chCoef[2][1], 2) / (2. * (variance + count*MODIFIER)));
					app[7] *= exp(-pow(rx[i][1] + chCoef[0][1] + chCoef[1][1] + chCoef[2][1], 2) / (2. * (variance + count*MODIFIER)));
				}
				else
				{
					app[0] *= exp(-pow(rx[i][1] - estimate[0][1] - estimate[1][1] - estimate[2][1], 2) / (2. * (variance + count*MODIFIER)));
					app[1] *= exp(-pow(rx[i][1] - estimate[0][1] - estimate[1][1] + estimate[2][1], 2) / (2. * (variance + count*MODIFIER)));
					app[2] *= exp(-pow(rx[i][1] - estimate[0][1] + estimate[1][1] - estimate[2][1], 2) / (2. * (variance + count*MODIFIER)));
					app[3] *= exp(-pow(rx[i][1] - estimate[0][1] + estimate[1][1] + estimate[2][1], 2) / (2. * (variance + count*MODIFIER)));
					app[4] *= exp(-pow(rx[i][1] + estimate[0][1] - estimate[1][1] - estimate[2][1], 2) / (2. * (variance + count*MODIFIER)));
					app[5] *= exp(-pow(rx[i][1] + estimate[0][1] - estimate[1][1] + estimate[2][1], 2) / (2. * (variance + count*MODIFIER)));
					app[6] *= exp(-pow(rx[i][1] + estimate[0][1] + estimate[1][1] - estimate[2][1], 2) / (2. * (variance + count*MODIFIER)));
					app[7] *= exp(-pow(rx[i][1] + estimate[0][1] + estimate[1][1] + estimate[2][1], 2) / (2. * (variance + count*MODIFIER)));
				}
			}
			if (MODIFIED_LLR)
			{
				for (int j = 0; j < NUM_LEVEL; j++)
				{
					if (app[j] == 0)
					{
						flag = false;
						break;
					}
				}
				if (!flag)
				{
					count++;
					continue;
				}
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
			if (flag) break;
		}
		//---------- a-posteriori LLRs ----------
		if (NUM_USER == 1)
		{
			if (app[0] <= NUMERIC_LIMIT) appLlr[0][i] = -LLR_LIMIT;
			else if (app[1] <= NUMERIC_LIMIT) appLlr[0][i] = LLR_LIMIT; 
			else appLlr[0][i] = log(app[0] / app[1]); 
		}
		else if (NUM_USER == 2)
		{
			//---------- user-A ----------
			if ((app[0] + app[1]) <= NUMERIC_LIMIT) appLlr[0][i] = -LLR_LIMIT;
			else if ((app[2] + app[3]) <= NUMERIC_LIMIT) appLlr[0][i] = LLR_LIMIT;
			else appLlr[0][i] = log((app[0] + app[1]) / (app[2] + app[3]));
			//---------- user-B ----------
			if ((app[0] + app[2]) <= NUMERIC_LIMIT) appLlr[1][i] = -LLR_LIMIT;
			else if ((app[1] + app[3]) <= NUMERIC_LIMIT) appLlr[1][i] = LLR_LIMIT;
			else appLlr[1][i] = log((app[0] + app[2]) / (app[1] + app[3]));
		}
		else if (NUM_USER == 3)
		{
			//---------- user-A ----------
			if ((app[0] + app[1] + app[2] + app[3]) <= NUMERIC_LIMIT) appLlr[0][i] = -LLR_LIMIT;
			else if ((app[4] + app[5] + app[6] + app[7]) <= NUMERIC_LIMIT) appLlr[0][i] = LLR_LIMIT;
			else appLlr[0][i] = log((app[0] + app[1] + app[2] + app[3]) / (app[4] + app[5] + app[6] + app[7]));
			//---------- user-B ----------
			if ((app[0] + app[1] + app[4] + app[5]) <= NUMERIC_LIMIT) appLlr[1][i] = -LLR_LIMIT;
			else if ((app[2] + app[3] + app[6] + app[7]) <= NUMERIC_LIMIT) appLlr[1][i] = LLR_LIMIT;
			else appLlr[1][i] = log((app[0] + app[1] + app[4] + app[5]) / (app[2] + app[3] + app[6] + app[7]));
			//---------- user-C ----------
			if ((app[0] + app[2] + app[4] + app[6]) <= NUMERIC_LIMIT) appLlr[2][i] = -LLR_LIMIT;
			else if ((app[1] + app[3] + app[5] + app[7]) <= NUMERIC_LIMIT) appLlr[2][i] = LLR_LIMIT;
			else appLlr[2][i] = log((app[0] + app[2] + app[4] + app[6]) / (app[1] + app[3] + app[5] + app[7]));
		}
	}
}