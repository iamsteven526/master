#define _USE_MATH_DEFINES
#include <iostream>
#include "parameters.h"
#include <cstring>
#include <math.h>
using namespace std;

void MLDT(double variance, double **chCoef, double ****supLevel, double ***postRx, double ***app, double **appLlr)
{
	for (int j = 0; j < NUM_USER; j++)
	{
		memset(appLlr[j], 0, sizeof(double)*BLOCK_LEN);
	}
	for(int m = 0; m < BLOCK_LEN/2; ++m){
		for (int i = 0; i < NUM_TX; i++)
		{
			//---------- a-posteriori prob. ----------
			if (NUM_USER == 1)
			{
				app[0][2*m + i][0] = exp(-pow(postRx[0][2*m + i][0] - supLevel[0][i][0][0], 2) / (2. * variance));
				app[0][2*m + i][1] = exp(-pow(postRx[0][2*m + i][0] + supLevel[0][i][0][0], 2) / (2. * variance));
			}
			else if (NUM_USER == 2)
			{
				for (int j = 0; j < NUM_USER; j++)
				{
					//cout << chCoef[0][0] << "  " << chCoef[0][1] << "  " << chCoef[1][0] << "  " << chCoef[1][1] << "  " << chCoef[2][0] << "  " << chCoef[2][1] << "  " << chCoef[3][0] << "  " << chCoef[3][1] << endl;
					//cout << postRx[j][2*m  + i][0] << "   " << supLevel[j][i][0][0] << "  " << supLevel[j][i][1][0] << "   " << supLevel[j][i][2][0] << endl;
					app[j][2*m  + i][0] = exp(-pow(postRx[j][2*m  + i][0] - supLevel[j][i][0][0] - supLevel[j][i][1][0] - supLevel[j][i][2][0], 2) / (2. * variance));
					app[j][2*m  + i][1] = exp(-pow(postRx[j][2*m  + i][0] - supLevel[j][i][0][0] - supLevel[j][i][1][0] + supLevel[j][i][2][0], 2) / (2. * variance));
					app[j][2*m  + i][2] = exp(-pow(postRx[j][2*m  + i][0] - supLevel[j][i][0][0] + supLevel[j][i][1][0] - supLevel[j][i][2][0], 2) / (2. * variance));
					app[j][2*m  + i][3] = exp(-pow(postRx[j][2*m  + i][0] - supLevel[j][i][0][0] + supLevel[j][i][1][0] + supLevel[j][i][2][0], 2) / (2. * variance));
					app[j][2*m  + i][4] = exp(-pow(postRx[j][2*m  + i][0] + supLevel[j][i][0][0] - supLevel[j][i][1][0] - supLevel[j][i][2][0], 2) / (2. * variance));
					app[j][2*m  + i][5] = exp(-pow(postRx[j][2*m  + i][0] + supLevel[j][i][0][0] - supLevel[j][i][1][0] + supLevel[j][i][2][0], 2) / (2. * variance));
					app[j][2*m  + i][6] = exp(-pow(postRx[j][2*m  + i][0] + supLevel[j][i][0][0] + supLevel[j][i][1][0] - supLevel[j][i][2][0], 2) / (2. * variance));
					app[j][2*m  + i][7] = exp(-pow(postRx[j][2*m  + i][0] + supLevel[j][i][0][0] + supLevel[j][i][1][0] + supLevel[j][i][2][0], 2) / (2. * variance));
                    //cout << app[j][2*m  + i][0] << "  " << app[j][2*m  + i][1] << "  " << app[j][2*m  + i][2] << "  " << app[j][2*m  + i][3] << "   " << app[j][2*m  + i][4] << "  " << app[j][2*m  + i][5] << "  " << app[j][2*m  + i][6] << "  " << app[j][2*m  + i][7] << "  "<< endl ;
				}
			}
			//---------- normalization ----------
			for (int j = 0; j < NUM_USER; j++)
			{
				double temp = 0;
				for (int k = 0; k < NUM_LEVEL/2; k++)
				{
					if (app[j][2*m + i][k] < NUMERIC_LIMIT) app[j][2*m  + i][k] = NUMERIC_LIMIT;
					temp += app[j][2*m  + i][k];
				}
				for (int k = 0; k < NUM_LEVEL/2; k++)
				{
					app[j][2*m  + i][k] /= temp;
					if (app[j][2*m  + i][k] < NUMERIC_LIMIT) app[j][2*m  + i][k] = NUMERIC_LIMIT;
					//cout << j << "  " << m << "   " << i << "  " << app[j][2*m + i][k] << endl;
				}
			}
			//---------- a-posteriori prob. ----------
			if (NUM_USER == 1)
			{
				app[0][2*m  + i][0] *= exp(-pow(postRx[0][2*m  + i][1] - supLevel[0][i][0][1], 2) / (2. * variance));
				app[0][2*m  + i][1] *= exp(-pow(postRx[0][2*m  + i][1] + supLevel[0][i][0][1], 2) / (2. * variance));
			}
			else if (NUM_USER == 2)
			{
				for (int j = 0; j < NUM_USER; j++)
				{
					app[j][2*m  + i][0] *= exp(-pow(postRx[j][2*m  + i][1] - supLevel[j][i][0][1] - supLevel[j][i][1][1] - supLevel[j][i][2][1], 2) / (2. * variance));
					app[j][2*m  + i][1] *= exp(-pow(postRx[j][2*m  + i][1] - supLevel[j][i][0][1] - supLevel[j][i][1][1] + supLevel[j][i][2][1], 2) / (2. * variance));
					app[j][2*m  + i][2] *= exp(-pow(postRx[j][2*m  + i][1] - supLevel[j][i][0][1] + supLevel[j][i][1][1] - supLevel[j][i][2][1], 2) / (2. * variance));
					app[j][2*m  + i][3] *= exp(-pow(postRx[j][2*m  + i][1] - supLevel[j][i][0][1] + supLevel[j][i][1][1] + supLevel[j][i][2][1], 2) / (2. * variance));
					app[j][2*m  + i][4] *= exp(-pow(postRx[j][2*m  + i][1] + supLevel[j][i][0][1] - supLevel[j][i][1][1] - supLevel[j][i][2][1], 2) / (2. * variance));
					app[j][2*m  + i][5] *= exp(-pow(postRx[j][2*m  + i][1] + supLevel[j][i][0][1] - supLevel[j][i][1][1] + supLevel[j][i][2][1], 2) / (2. * variance));
					app[j][2*m  + i][6] *= exp(-pow(postRx[j][2*m  + i][1] + supLevel[j][i][0][1] + supLevel[j][i][1][1] - supLevel[j][i][2][1], 2) / (2. * variance));
					app[j][2*m  + i][7] *= exp(-pow(postRx[j][2*m  + i][1] + supLevel[j][i][0][1] + supLevel[j][i][1][1] + supLevel[j][i][2][1], 2) / (2. * variance));
				}
			}
			//---------- normalization ----------
			for (int j = 0; j < NUM_USER; j++)
			{
				double temp = 0;
				for (int k = 0; k < NUM_LEVEL/2; k++)
				{
					if (app[j][2*m  + i][k] < NUMERIC_LIMIT) app[j][2*m  + i][k] = NUMERIC_LIMIT;
					temp += app[j][2*m  + i][k];
				}
				for (int k = 0; k < NUM_LEVEL/2; k++)
				{
					app[j][2*m  + i][k] /= temp;
					if (app[j][2*m  + i][k] < NUMERIC_LIMIT) app[j][2*m  + i][k] = NUMERIC_LIMIT;
					//cout << app[j][2*m + i][k] << endl;
				}
			}
		}
		//---------- a-posteriori LLRs ----------
		if (NUM_USER == 1)
		{
			for (int i = 0; i < NUM_TX; i++)
			{
				if (app[0][2*m  + i][0] <= NUMERIC_LIMIT) appLlr[0][2*m  + i] = -LLR_LIMIT;
				else if (app[0][2*m  + i][1] <= NUMERIC_LIMIT) appLlr[0][2*m  + i] = LLR_LIMIT;
				else appLlr[0][2*m  + i] = log(app[0][2*m  + i][0] / app[0][2*m  + i][1]);
			}
		}
		else if (NUM_USER == 2)
		{
			for (int j = 0; j < NUM_USER; j++)
			{
				for (int i = 0; i < NUM_TX; i++)
				{
					//---------- desired signal ----------
					if ((app[j][2*m  + i][0] + app[j][2*m  + i][1] + app[j][2*m  + i][2] + app[j][2*m  + i][3]) <= NUMERIC_LIMIT) appLlr[j][2*m  + i] += -LLR_LIMIT;
					else if ((app[j][2*m  + i][4] + app[j][2*m  + i][5] + app[j][2*m  + i][6] + app[j][2*m  + i][7]) <= NUMERIC_LIMIT) appLlr[j][2*m  + i] += LLR_LIMIT;
					else appLlr[j][2*m  + i] += log((app[j][2*m  + i][0] + app[j][2*m  + i][1] + app[j][2*m  + i][2] + app[j][2*m  + i][3]) / (app[j][2*m  + i][4] + app[j][2*m  + i][5] + app[j][2*m  + i][6] + app[j][2*m  + i][7]));
					//---------- interferences ----------
					if ((app[j][2*m  + i][0] + app[j][2*m  + i][1] + app[j][2*m  + i][4] + app[j][2*m  + i][5]) <= NUMERIC_LIMIT) appLlr[(j + 1) % NUM_USER][2*m] += -LLR_LIMIT;
					else if ((app[j][2*m  + i][2] + app[j][2*m  + i][3] + app[j][2*m  + i][6] + app[j][2*m  + i][7]) <= NUMERIC_LIMIT) appLlr[(j + 1) % NUM_USER][2*m] += LLR_LIMIT;
					else appLlr[(j + 1) % NUM_USER][2*m] += log((app[j][2*m  + i][0] + app[j][2*m  + i][1] + app[j][2*m  + i][4] + app[j][2*m  + i][5]) / (app[j][2*m  + i][2] + app[j][2*m  + i][3] + app[j][2*m  + i][6] + app[j][2*m  + i][7]));
					if ((app[j][2*m  + i][0] + app[j][2*m  + i][2] + app[j][2*m  + i][4] + app[j][2*m  + i][6]) <= NUMERIC_LIMIT) appLlr[(j + 1) % NUM_USER][2*m +1] += -LLR_LIMIT;
					else if ((app[j][2*m  + i][1] + app[j][2*m  + i][3] + app[j][2*m  + i][5] + app[j][2*m  + i][7]) <= NUMERIC_LIMIT) appLlr[(j + 1) % NUM_USER][2*m + 1] += LLR_LIMIT;
					else appLlr[(j + 1) % NUM_USER][2*m + 1] += log((app[j][2*m  + i][0] + app[j][2*m  + i][2] + app[j][2*m  + i][4] + app[j][2*m  + i][6]) / (app[j][2*m  + i][1] + app[j][2*m  + i][3] + app[j][2*m  + i][5] + app[j][2*m  + i][7]));
				    //cout << m << "   " << appLlr[j][2*m + i] << endl;
				}
			}
		}
	}
}