#define _USE_MATH_DEFINES
#include <iostream>
#include "parameters.h"
#include <cstring>
#include <math.h>
using namespace std;

void MLDT(double variance, double ***chCoef, double ****supLevel, double ***postRx, double ***app, double **appLlr)
{
	for (int i = 0; i < NUM_TX; i++)
	{
		//---------- a-posteriori prob. ----------
		if (NUM_USER == 1)
		{
			app[0][i][0] = exp(-pow(postRx[0][i][0] - supLevel[0][i][0][0], 2) / (2. * variance));
			app[0][i][1] = exp(-pow(postRx[0][i][0] + supLevel[0][i][0][0], 2) / (2. * variance));
		}
		else if (NUM_USER == 2)
		{
			for (int j = 0; j < NUM_USER; j++)
			{
				app[j][i][0] = exp(-pow(postRx[j][i][0] - supLevel[j][i][0][0] - supLevel[j][i][1][0] - supLevel[j][i][2][0], 2) / (2. * variance));
				app[j][i][1] = exp(-pow(postRx[j][i][0] - supLevel[j][i][0][0] - supLevel[j][i][1][0] + supLevel[j][i][2][0], 2) / (2. * variance));
				app[j][i][2] = exp(-pow(postRx[j][i][0] - supLevel[j][i][0][0] + supLevel[j][i][1][0] - supLevel[j][i][2][0], 2) / (2. * variance));
				app[j][i][3] = exp(-pow(postRx[j][i][0] - supLevel[j][i][0][0] + supLevel[j][i][1][0] + supLevel[j][i][2][0], 2) / (2. * variance));
				app[j][i][4] = exp(-pow(postRx[j][i][0] + supLevel[j][i][0][0] - supLevel[j][i][1][0] - supLevel[j][i][2][0], 2) / (2. * variance));
				app[j][i][5] = exp(-pow(postRx[j][i][0] + supLevel[j][i][0][0] - supLevel[j][i][1][0] + supLevel[j][i][2][0], 2) / (2. * variance));
				app[j][i][6] = exp(-pow(postRx[j][i][0] + supLevel[j][i][0][0] + supLevel[j][i][1][0] - supLevel[j][i][2][0], 2) / (2. * variance));
				app[j][i][7] = exp(-pow(postRx[j][i][0] + supLevel[j][i][0][0] + supLevel[j][i][1][0] + supLevel[j][i][2][0], 2) / (2. * variance));
			}
		}
		else if (NUM_USER == 3)
		{
			for (int j = 0; j < NUM_USER; j++)
			{
				for (int k = 0; k < 32; ++k){
					app[j][i][k] = exp(-pow(postRx[j][i][0] - pow(-1,k/16)*supLevel[j][i][0][0] - pow(-1,k/8)*supLevel[j][i][1][0] - pow(-1,k/4)*supLevel[j][i][2][0]  - pow(-1,k/2)*supLevel[j][i][3][0] - pow(-1,k)*supLevel[j][i][4][0], 2) / (2. * variance));
				}
			}
		}
		//---------- normalization ----------
		for (int j = 0; j < NUM_USER; j++)
		{
			double temp = 0;
			for (int k = 0; k < NUM_LEVEL; k++)
			{
				if (app[j][i][k] < NUMERIC_LIMIT) app[j][i][k] = NUMERIC_LIMIT;
				temp += app[j][i][k];
			}
			for (int k = 0; k < NUM_LEVEL; k++)
			{
				app[j][i][k] /= temp;
				if (app[j][i][k] < NUMERIC_LIMIT) app[j][i][k] = NUMERIC_LIMIT;
			}
		}
		//---------- a-posteriori prob. ----------
		if (NUM_USER == 1)
		{
			app[0][i][0] *= exp(-pow(postRx[0][i][1] - supLevel[0][i][0][1], 2) / (2. * variance));
			app[0][i][1] *= exp(-pow(postRx[0][i][1] + supLevel[0][i][0][1], 2) / (2. * variance));
		}
		else if (NUM_USER == 2)
		{
			for (int j = 0; j < NUM_USER; j++)
			{
				app[j][i][0] *= exp(-pow(postRx[j][i][1] - supLevel[j][i][0][1] - supLevel[j][i][1][1] - supLevel[j][i][2][1], 2) / (2. * variance));
				app[j][i][1] *= exp(-pow(postRx[j][i][1] - supLevel[j][i][0][1] - supLevel[j][i][1][1] + supLevel[j][i][2][1], 2) / (2. * variance));
				app[j][i][2] *= exp(-pow(postRx[j][i][1] - supLevel[j][i][0][1] + supLevel[j][i][1][1] - supLevel[j][i][2][1], 2) / (2. * variance));
				app[j][i][3] *= exp(-pow(postRx[j][i][1] - supLevel[j][i][0][1] + supLevel[j][i][1][1] + supLevel[j][i][2][1], 2) / (2. * variance));
				app[j][i][4] *= exp(-pow(postRx[j][i][1] + supLevel[j][i][0][1] - supLevel[j][i][1][1] - supLevel[j][i][2][1], 2) / (2. * variance));
				app[j][i][5] *= exp(-pow(postRx[j][i][1] + supLevel[j][i][0][1] - supLevel[j][i][1][1] + supLevel[j][i][2][1], 2) / (2. * variance));
				app[j][i][6] *= exp(-pow(postRx[j][i][1] + supLevel[j][i][0][1] + supLevel[j][i][1][1] - supLevel[j][i][2][1], 2) / (2. * variance));
				app[j][i][7] *= exp(-pow(postRx[j][i][1] + supLevel[j][i][0][1] + supLevel[j][i][1][1] + supLevel[j][i][2][1], 2) / (2. * variance));
			}
		}
		else if (NUM_USER == 3)
		{
			for (int j = 0; j < NUM_USER; j++)
			{
				for (int k = 0; k < 32; ++k){
					app[j][i][k] *= exp(-pow(postRx[j][i][1] - pow(-1,k/16)*supLevel[j][i][0][1] - pow(-1,k/8)*supLevel[j][i][1][1] - pow(-1,k/4)*supLevel[j][i][2][1]  - pow(-1,k/2)*supLevel[j][i][3][1] - pow(-1,k)*supLevel[j][i][4][1], 2) / (2. * variance));
				}
			}
		}
		//---------- normalization ----------
		for (int j = 0; j < NUM_USER; j++)
		{
			double temp = 0;
			for (int k = 0; k < NUM_LEVEL; k++)
			{
				if (app[j][i][k] < NUMERIC_LIMIT) app[j][i][k] = NUMERIC_LIMIT;
				temp += app[j][i][k];
			}
			for (int k = 0; k < NUM_LEVEL; k++)
			{
				app[j][i][k] /= temp;
				if (app[j][i][k] < NUMERIC_LIMIT) app[j][i][k] = NUMERIC_LIMIT;
			}
		}
	}
	//---------- a-posteriori LLRs ----------
	if (NUM_USER == 1)
	{
		for (int i = 0; i < NUM_TX; i++)
		{
			if (app[0][i][0] <= NUMERIC_LIMIT) appLlr[0][i] = -LLR_LIMIT;
			else if (app[0][i][1] <= NUMERIC_LIMIT) appLlr[0][i] = LLR_LIMIT;
			else appLlr[0][i] = log(app[0][i][0] / app[0][i][1]);
		}
	}
	else if (NUM_USER == 2)
	{
		for (int j = 0; j < NUM_USER; j++)
		{
			memset(appLlr[j], 0, sizeof(double)*NUM_TX);
		}
		for (int j = 0; j < NUM_USER; j++)
		{
			for (int i = 0; i < NUM_TX; i++)
			{
				//---------- desired signal ----------
				if ((app[j][i][0] + app[j][i][1] + app[j][i][2] + app[j][i][3]) <= NUMERIC_LIMIT) appLlr[j][i] += -2*LLR_LIMIT;
				else if ((app[j][i][4] + app[j][i][5] + app[j][i][6] + app[j][i][7]) <= NUMERIC_LIMIT) appLlr[j][i] += 2*LLR_LIMIT;
				else appLlr[j][i] += log((app[j][i][0] + app[j][i][1] + app[j][i][2] + app[j][i][3]) / (app[j][i][4] + app[j][i][5] + app[j][i][6] + app[j][i][7]));
				//---------- interferences ----------
				if ((app[j][i][0] + app[j][i][1] + app[j][i][4] + app[j][i][5]) <= NUMERIC_LIMIT) appLlr[(j + 1) % NUM_USER][0] += -LLR_LIMIT;
				else if ((app[j][i][2] + app[j][i][3] + app[j][i][6] + app[j][i][7]) <= NUMERIC_LIMIT) appLlr[(j + 1) % NUM_USER][0] += LLR_LIMIT;
				else appLlr[(j + 1) % NUM_USER][0] += log((app[j][i][0] + app[j][i][1] + app[j][i][4] + app[j][i][5]) / (app[j][i][2] + app[j][i][3] + app[j][i][6] + app[j][i][7]));
				if ((app[j][i][0] + app[j][i][2] + app[j][i][4] + app[j][i][6]) <= NUMERIC_LIMIT) appLlr[(j + 1) % NUM_USER][1] += -LLR_LIMIT;
				else if ((app[j][i][1] + app[j][i][3] + app[j][i][5] + app[j][i][7]) <= NUMERIC_LIMIT) appLlr[(j + 1) % NUM_USER][1] += LLR_LIMIT;
				else appLlr[(j + 1) % NUM_USER][1] += log((app[j][i][0] + app[j][i][2] + app[j][i][4] + app[j][i][6]) / (app[j][i][1] + app[j][i][3] + app[j][i][5] + app[j][i][7]));
			}
		}
	}
	else if (NUM_USER == 3)
	{
		for (int j = 0; j < NUM_USER; j++)
		{
			memset(appLlr[j], 0, sizeof(double)*NUM_TX);
		}	
		for (int j = 0; j < NUM_USER; j++)
		{
		    for (int i = 0; i < NUM_TX; i++)
			{
				//TODO:::::
				double prob0, prob1;
				//---------- desired signal ---------- (1)
				prob0 = app[j][i][0] + app[j][i][1] + app[j][i][2] + app[j][i][3] + app[j][i][4] + app[j][i][5] + app[j][i][6] + app[j][i][7] + app[j][i][8] + app[j][i][9] + app[j][i][10] + app[j][i][11] + app[j][i][12] + app[j][i][13] + app[j][i][14] + app[j][i][15];
				prob1 = app[j][i][16] + app[j][i][17] + app[j][i][18] + app[j][i][19] + app[j][i][20] + app[j][i][21] + app[j][i][22] + app[j][i][23] + app[j][i][24] + app[j][i][25] + app[j][i][26] + app[j][i][27] + app[j][i][28] + app[j][i][29] + app[j][i][30] + app[j][i][31];
				if (prob0 <= NUMERIC_LIMIT) appLlr[j][i] += -2*LLR_LIMIT;
				else if (prob1 <= NUMERIC_LIMIT) appLlr[j][i] += 2*LLR_LIMIT;
				else appLlr[j][i] += log(prob0 / prob1);
				//---------- interferences ----------  (4)
				/*
				prob0 = app[j][i][0] + app[j][i][1] + app[j][i][2] + app[j][i][3] + app[j][i][4] + app[j][i][5] + app[j][i][6] + app[j][i][7] + app[j][i][16] + app[j][i][17] + app[j][i][18] + app[j][i][19] + app[j][i][20] + app[j][i][21] + app[j][i][22] + app[j][i][23];
				prob1 = app[j][i][8] + app[j][i][9] + app[j][i][10] + app[j][i][11] + app[j][i][12] + app[j][i][13] + app[j][i][14] + app[j][i][15] + app[j][i][24] + app[j][i][25] + app[j][i][26] + app[j][i][27] + app[j][i][28] + app[j][i][29] + app[j][i][30] + app[j][i][31];
				if (prob0 <= NUMERIC_LIMIT) appLlr[j][0] += -2*LLR_LIMIT;
				else if (prob1 <= NUMERIC_LIMIT) appLlr[j][0] += 2*LLR_LIMIT;
				else appLlr[j][0] += log(prob0 / prob1);
				*/






                /*
				if ((app[j][i][0] + app[j][i][1] + app[j][i][4] + app[j][i][5]) <= NUMERIC_LIMIT) appLlr[(j + 1) % NUM_USER][0] += -LLR_LIMIT;
				else if ((app[j][i][2] + app[j][i][3] + app[j][i][6] + app[j][i][7]) <= NUMERIC_LIMIT) appLlr[(j + 1) % NUM_USER][0] += LLR_LIMIT;
				else appLlr[(j + 1) % NUM_USER][0] += log((app[j][i][0] + app[j][i][1] + app[j][i][4] + app[j][i][5]) / (app[j][i][2] + app[j][i][3] + app[j][i][6] + app[j][i][7]));
				if ((app[j][i][0] + app[j][i][2] + app[j][i][4] + app[j][i][6]) <= NUMERIC_LIMIT) appLlr[(j + 1) % NUM_USER][1] += -LLR_LIMIT;
				else if ((app[j][i][1] + app[j][i][3] + app[j][i][5] + app[j][i][7]) <= NUMERIC_LIMIT) appLlr[(j + 1) % NUM_USER][1] += LLR_LIMIT;
				else appLlr[(j + 1) % NUM_USER][1] += log((app[j][i][0] + app[j][i][2] + app[j][i][4] + app[j][i][6]) / (app[j][i][1] + app[j][i][3] + app[j][i][5] + app[j][i][7]));
			    */
			}
		}	
	}
}