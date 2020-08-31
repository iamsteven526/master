#define _USE_MATH_DEFINES
#include <iostream>
#include "parameters.h"
#include <cmath>
#include <cstring>
using namespace std;

void MLDT(LDPC &ldpc, double variance, double ****H, double ***postRx, double **app, double **appLlr, double **refLlr, double ****estimate)
{
	//---- Different P/S for Differential encoder and NBC encoder
	
	if (CH_CODING_TYPE)
	{
		for (int i = 0, m = 0; i < FFT_SEGMENT; i++)
		{
			for (int j = 0; j < FFT_POINT; j++)
			{
				//---------- APPs ----------
				if (NUM_USER == 1)
				{
					if (CE_SCHEME == 1)
					{
						app[m][0] = exp(-pow(postRx[i][0][j] - H[0][i][0][j], 2) / (2. * variance));
						app[m][1] = exp(-pow(postRx[i][0][j] + H[0][i][0][j], 2) / (2. * variance));
					}
					else
					{
						app[m][0] = exp(-pow(postRx[i][0][j] - estimate[0][i * SLIDING][j][0], 2) / (2. * variance));
						app[m][1] = exp(-pow(postRx[i][0][j] + estimate[0][i * SLIDING][j][0], 2) / (2. * variance));
					}
				}
				else if (NUM_USER == 2)
				{
					if (CE_SCHEME == 1)
					{
						app[m][0] = exp(-pow(postRx[i][0][j] - H[0][i][0][j] - H[1][i][0][j], 2) / (2. * variance));
						app[m][1] = exp(-pow(postRx[i][0][j] + H[0][i][0][j] - H[1][i][0][j], 2) / (2. * variance));
						app[m][2] = exp(-pow(postRx[i][0][j] - H[0][i][0][j] + H[1][i][0][j], 2) / (2. * variance));
						app[m][3] = exp(-pow(postRx[i][0][j] + H[0][i][0][j] + H[1][i][0][j], 2) / (2. * variance));
					}
					else
					{
						app[m][0] = exp(-pow(postRx[i][0][j] - estimate[0][i * SLIDING][j][0] - estimate[1][i * SLIDING][j][0], 2) / (2. * variance));
						app[m][1] = exp(-pow(postRx[i][0][j] + estimate[0][i * SLIDING][j][0] - estimate[1][i * SLIDING][j][0], 2) / (2. * variance));
						app[m][2] = exp(-pow(postRx[i][0][j] - estimate[0][i * SLIDING][j][0] + estimate[1][i * SLIDING][j][0], 2) / (2. * variance));
						app[m][3] = exp(-pow(postRx[i][0][j] + estimate[0][i * SLIDING][j][0] + estimate[1][i * SLIDING][j][0], 2) / (2. * variance));
					}
				}
				else if (NUM_USER == 3)
				{
					if (CE_SCHEME == 1)
					{
						app[m][0] = exp(-pow(postRx[i][0][j] - H[0][i][0][j] - H[1][i][0][j] - H[2][i][0][j], 2) / (2. * variance));
						app[m][1] = exp(-pow(postRx[i][0][j] - H[0][i][0][j] - H[1][i][0][j] + H[2][i][0][j], 2) / (2. * variance));
						app[m][2] = exp(-pow(postRx[i][0][j] - H[0][i][0][j] + H[1][i][0][j] - H[2][i][0][j], 2) / (2. * variance));
						app[m][3] = exp(-pow(postRx[i][0][j] - H[0][i][0][j] + H[1][i][0][j] + H[2][i][0][j], 2) / (2. * variance));
						app[m][4] = exp(-pow(postRx[i][0][j] + H[0][i][0][j] - H[1][i][0][j] - H[2][i][0][j], 2) / (2. * variance));
						app[m][5] = exp(-pow(postRx[i][0][j] + H[0][i][0][j] - H[1][i][0][j] + H[2][i][0][j], 2) / (2. * variance));
						app[m][6] = exp(-pow(postRx[i][0][j] + H[0][i][0][j] + H[1][i][0][j] - H[2][i][0][j], 2) / (2. * variance));
						app[m][7] = exp(-pow(postRx[i][0][j] + H[0][i][0][j] + H[1][i][0][j] + H[2][i][0][j], 2) / (2. * variance));
					}
					else
					{
						app[m][0] = exp(-pow(postRx[i][0][j] - estimate[0][i * SLIDING][j][0] - estimate[1][i * SLIDING][j][0] - estimate[2][i * SLIDING][j][0], 2) / (2. * variance));
						app[m][1] = exp(-pow(postRx[i][0][j] - estimate[0][i * SLIDING][j][0] - estimate[1][i * SLIDING][j][0] + estimate[2][i * SLIDING][j][0], 2) / (2. * variance));
						app[m][2] = exp(-pow(postRx[i][0][j] - estimate[0][i * SLIDING][j][0] + estimate[1][i * SLIDING][j][0] - estimate[2][i * SLIDING][j][0], 2) / (2. * variance));
						app[m][3] = exp(-pow(postRx[i][0][j] - estimate[0][i * SLIDING][j][0] + estimate[1][i * SLIDING][j][0] + estimate[2][i * SLIDING][j][0], 2) / (2. * variance));
						app[m][4] = exp(-pow(postRx[i][0][j] + estimate[0][i * SLIDING][j][0] - estimate[1][i * SLIDING][j][0] - estimate[2][i * SLIDING][j][0], 2) / (2. * variance));
						app[m][5] = exp(-pow(postRx[i][0][j] + estimate[0][i * SLIDING][j][0] - estimate[1][i * SLIDING][j][0] + estimate[2][i * SLIDING][j][0], 2) / (2. * variance));
						app[m][6] = exp(-pow(postRx[i][0][j] + estimate[0][i * SLIDING][j][0] + estimate[1][i * SLIDING][j][0] - estimate[2][i * SLIDING][j][0], 2) / (2. * variance));
						app[m][7] = exp(-pow(postRx[i][0][j] + estimate[0][i * SLIDING][j][0] + estimate[1][i * SLIDING][j][0] + estimate[2][i * SLIDING][j][0], 2) / (2. * variance));
					}
				}
				else
				{
					printf("\nPARAMETER SETTING IS WRONG\n");
					system("pause");
				}
				//---------- normalization ----------
				double temp = 0;
				for (int k = 0; k < NUM_LEVEL; k++)
				{
					if (app[m][k] < NUMERIC_LIMIT) app[m][k] = NUMERIC_LIMIT;
					temp += app[m][k];
				}
				for (int k = 0; k < NUM_LEVEL; k++)
				{
					app[m][k] /= temp;
					if (app[m][k] < NUMERIC_LIMIT) app[m][k] = NUMERIC_LIMIT;
				}
				//---------- APPs ----------
				if (NUM_USER == 1)
				{
					if (CE_SCHEME == 1)
					{
						app[m][0] *= exp(-pow(postRx[i][1][j] - H[0][i][1][j], 2) / (2. * variance));
						app[m][1] *= exp(-pow(postRx[i][1][j] + H[0][i][1][j], 2) / (2. * variance));
					}
					else
					{
						app[m][0] *= exp(-pow(postRx[i][1][j] - estimate[0][i * SLIDING][j][1], 2) / (2. * variance));
						app[m][1] *= exp(-pow(postRx[i][1][j] + estimate[0][i * SLIDING][j][1], 2) / (2. * variance));
					}
				}
				else if (NUM_USER == 2)
				{
					if (CE_SCHEME == 1)
					{
						app[m][0] *= exp(-pow(postRx[i][1][j] - H[0][i][1][j] - H[1][i][1][j], 2) / (2. * variance));
						app[m][1] *= exp(-pow(postRx[i][1][j] + H[0][i][1][j] - H[1][i][1][j], 2) / (2. * variance));
						app[m][2] *= exp(-pow(postRx[i][1][j] - H[0][i][1][j] + H[1][i][1][j], 2) / (2. * variance));
						app[m][3] *= exp(-pow(postRx[i][1][j] + H[0][i][1][j] + H[1][i][1][j], 2) / (2. * variance));
					}
					else
					{
						app[m][0] *= exp(-pow(postRx[i][1][j] - estimate[0][i * SLIDING][j][1] - estimate[1][i * SLIDING][j][1], 2) / (2. * variance));
						app[m][1] *= exp(-pow(postRx[i][1][j] + estimate[0][i * SLIDING][j][1] - estimate[1][i * SLIDING][j][1], 2) / (2. * variance));
						app[m][2] *= exp(-pow(postRx[i][1][j] - estimate[0][i * SLIDING][j][1] + estimate[1][i * SLIDING][j][1], 2) / (2. * variance));
						app[m][3] *= exp(-pow(postRx[i][1][j] + estimate[0][i * SLIDING][j][1] + estimate[1][i * SLIDING][j][1], 2) / (2. * variance));
					}
				}
				else if (NUM_USER == 3)
				{
					if (CE_SCHEME == 1)
					{
						app[m][0] *= exp(-pow(postRx[i][1][j] - H[0][i][1][j] - H[1][i][1][j] - H[2][i][1][j], 2) / (2. * variance));
						app[m][1] *= exp(-pow(postRx[i][1][j] - H[0][i][1][j] - H[1][i][1][j] + H[2][i][1][j], 2) / (2. * variance));
						app[m][2] *= exp(-pow(postRx[i][1][j] - H[0][i][1][j] + H[1][i][1][j] - H[2][i][1][j], 2) / (2. * variance));
						app[m][3] *= exp(-pow(postRx[i][1][j] - H[0][i][1][j] + H[1][i][1][j] + H[2][i][1][j], 2) / (2. * variance));
						app[m][4] *= exp(-pow(postRx[i][1][j] + H[0][i][1][j] - H[1][i][1][j] - H[2][i][1][j], 2) / (2. * variance));
						app[m][5] *= exp(-pow(postRx[i][1][j] + H[0][i][1][j] - H[1][i][1][j] + H[2][i][1][j], 2) / (2. * variance));
						app[m][6] *= exp(-pow(postRx[i][1][j] + H[0][i][1][j] + H[1][i][1][j] - H[2][i][1][j], 2) / (2. * variance));
						app[m][7] *= exp(-pow(postRx[i][1][j] + H[0][i][1][j] + H[1][i][1][j] + H[2][i][1][j], 2) / (2. * variance));
					}
					else
					{
						app[m][0] *= exp(-pow(postRx[i][1][j] - estimate[0][i * SLIDING][j][1] - estimate[1][i * SLIDING][j][1] - estimate[2][i * SLIDING][j][1], 2) / (2. * variance));
						app[m][1] *= exp(-pow(postRx[i][1][j] - estimate[0][i * SLIDING][j][1] - estimate[1][i * SLIDING][j][1] + estimate[2][i * SLIDING][j][1], 2) / (2. * variance));
						app[m][2] *= exp(-pow(postRx[i][1][j] - estimate[0][i * SLIDING][j][1] + estimate[1][i * SLIDING][j][1] - estimate[2][i * SLIDING][j][1], 2) / (2. * variance));
						app[m][3] *= exp(-pow(postRx[i][1][j] - estimate[0][i * SLIDING][j][1] + estimate[1][i * SLIDING][j][1] + estimate[2][i * SLIDING][j][1], 2) / (2. * variance));
						app[m][4] *= exp(-pow(postRx[i][1][j] + estimate[0][i * SLIDING][j][1] - estimate[1][i * SLIDING][j][1] - estimate[2][i * SLIDING][j][1], 2) / (2. * variance));
						app[m][5] *= exp(-pow(postRx[i][1][j] + estimate[0][i * SLIDING][j][1] - estimate[1][i * SLIDING][j][1] + estimate[2][i * SLIDING][j][1], 2) / (2. * variance));
						app[m][6] *= exp(-pow(postRx[i][1][j] + estimate[0][i * SLIDING][j][1] + estimate[1][i * SLIDING][j][1] - estimate[2][i * SLIDING][j][1], 2) / (2. * variance));
						app[m][7] *= exp(-pow(postRx[i][1][j] + estimate[0][i * SLIDING][j][1] + estimate[1][i * SLIDING][j][1] + estimate[2][i * SLIDING][j][1], 2) / (2. * variance));
					}
				}
				//---------- normalization ----------
				temp = 0;
				for (int k = 0; k < NUM_LEVEL; k++)
				{
					if (app[m][k] < NUMERIC_LIMIT) app[m][k] = NUMERIC_LIMIT;
					temp += app[m][k];
				}
				for (int k = 0; k < NUM_LEVEL; k++)
				{
					app[m][k] /= temp;
					if (app[m][k] < NUMERIC_LIMIT) app[m][k] = NUMERIC_LIMIT;
				}
				m++;
			}
		}
	}
	else
	{
		for (int j = 0, m = 0; j < FFT_POINT; j++)
		{
			for (int i = 0; i < FFT_SEGMENT; i++)
			{
				//---------- APPs ----------
				if (NUM_USER == 1)
				{
					if (CE_SCHEME == 1)
					{
						app[m][0] = exp(-pow(postRx[i][0][j] - H[0][i][0][j], 2) / (2. * variance));
						app[m][1] = exp(-pow(postRx[i][0][j] + H[0][i][0][j], 2) / (2. * variance));
					}
					else
					{
						app[m][0] = exp(-pow(postRx[i][0][j] - estimate[0][i * SLIDING][j][0], 2) / (2. * variance));
						app[m][1] = exp(-pow(postRx[i][0][j] + estimate[0][i * SLIDING][j][0], 2) / (2. * variance));
					}
				}
				else if (NUM_USER == 2)
				{

					if (CE_SCHEME == 1)
					{
						app[m][0] = exp(-pow(postRx[i][0][j] - H[0][i][0][j] - H[1][i][0][j], 2) / (2. * variance));
						app[m][1] = exp(-pow(postRx[i][0][j] - H[0][i][0][j] + H[1][i][0][j], 2) / (2. * variance));
						app[m][2] = exp(-pow(postRx[i][0][j] + H[0][i][0][j] - H[1][i][0][j], 2) / (2. * variance));
						app[m][3] = exp(-pow(postRx[i][0][j] + H[0][i][0][j] + H[1][i][0][j], 2) / (2. * variance));
					}
					else
					{

						app[m][0] = exp(-pow(postRx[i][0][j] - estimate[0][i * SLIDING][j][0] - estimate[1][i * SLIDING][j][0], 2) / (2. * variance));
						app[m][1] = exp(-pow(postRx[i][0][j] - estimate[0][i * SLIDING][j][0] + estimate[1][i * SLIDING][j][0], 2) / (2. * variance));
						app[m][2] = exp(-pow(postRx[i][0][j] + estimate[0][i * SLIDING][j][0] - estimate[1][i * SLIDING][j][0], 2) / (2. * variance));
						app[m][3] = exp(-pow(postRx[i][0][j] + estimate[0][i * SLIDING][j][0] + estimate[1][i * SLIDING][j][0], 2) / (2. * variance));
					}
				}
				else if (NUM_USER == 3)
				{
					if (CE_SCHEME == 1)
					{
						app[m][0] = exp(-pow(postRx[i][0][j] - H[0][i][0][j] - H[1][i][0][j] - H[2][i][0][j], 2) / (2. * variance));
						app[m][1] = exp(-pow(postRx[i][0][j] - H[0][i][0][j] - H[1][i][0][j] + H[2][i][0][j], 2) / (2. * variance));
						app[m][2] = exp(-pow(postRx[i][0][j] - H[0][i][0][j] + H[1][i][0][j] - H[2][i][0][j], 2) / (2. * variance));
						app[m][3] = exp(-pow(postRx[i][0][j] - H[0][i][0][j] + H[1][i][0][j] + H[2][i][0][j], 2) / (2. * variance));
						app[m][4] = exp(-pow(postRx[i][0][j] + H[0][i][0][j] - H[1][i][0][j] - H[2][i][0][j], 2) / (2. * variance));
						app[m][5] = exp(-pow(postRx[i][0][j] + H[0][i][0][j] - H[1][i][0][j] + H[2][i][0][j], 2) / (2. * variance));
						app[m][6] = exp(-pow(postRx[i][0][j] + H[0][i][0][j] + H[1][i][0][j] - H[2][i][0][j], 2) / (2. * variance));
						app[m][7] = exp(-pow(postRx[i][0][j] + H[0][i][0][j] + H[1][i][0][j] + H[2][i][0][j], 2) / (2. * variance));
					}
					else
					{
						app[m][0] = exp(-pow(postRx[i][0][j] - estimate[0][i * SLIDING][j][0] - estimate[1][i * SLIDING][j][0] - estimate[2][i * SLIDING][j][0], 2) / (2. * variance));
						app[m][1] = exp(-pow(postRx[i][0][j] - estimate[0][i * SLIDING][j][0] - estimate[1][i * SLIDING][j][0] + estimate[2][i * SLIDING][j][0], 2) / (2. * variance));
						app[m][2] = exp(-pow(postRx[i][0][j] - estimate[0][i * SLIDING][j][0] + estimate[1][i * SLIDING][j][0] - estimate[2][i * SLIDING][j][0], 2) / (2. * variance));
						app[m][3] = exp(-pow(postRx[i][0][j] - estimate[0][i * SLIDING][j][0] + estimate[1][i * SLIDING][j][0] + estimate[2][i * SLIDING][j][0], 2) / (2. * variance));
						app[m][4] = exp(-pow(postRx[i][0][j] + estimate[0][i * SLIDING][j][0] - estimate[1][i * SLIDING][j][0] - estimate[2][i * SLIDING][j][0], 2) / (2. * variance));
						app[m][5] = exp(-pow(postRx[i][0][j] + estimate[0][i * SLIDING][j][0] - estimate[1][i * SLIDING][j][0] + estimate[2][i * SLIDING][j][0], 2) / (2. * variance));
						app[m][6] = exp(-pow(postRx[i][0][j] + estimate[0][i * SLIDING][j][0] + estimate[1][i * SLIDING][j][0] - estimate[2][i * SLIDING][j][0], 2) / (2. * variance));
						app[m][7] = exp(-pow(postRx[i][0][j] + estimate[0][i * SLIDING][j][0] + estimate[1][i * SLIDING][j][0] + estimate[2][i * SLIDING][j][0], 2) / (2. * variance));
					}
				}
				else
				{
					printf("\nPARAMETER SETTING IS WRONG\n");
					system("pause");
				}
				//---------- normalization ----------
				double temp = 0;
				for (int k = 0; k < NUM_LEVEL; k++)
				{
					if (app[m][k] < NUMERIC_LIMIT) app[m][k] = NUMERIC_LIMIT;
					temp += app[m][k];
				}
				for (int k = 0; k < NUM_LEVEL; k++)
				{
					app[m][k] /= temp;
					if (app[m][k] < NUMERIC_LIMIT) app[m][k] = NUMERIC_LIMIT;
				}
				//---------- APPs ----------
				if (NUM_USER == 1)
				{
					if (CE_SCHEME == 1)
					{
						app[m][0] *= exp(-pow(postRx[i][1][j] - H[0][i][1][j], 2) / (2. * variance));
						app[m][1] *= exp(-pow(postRx[i][1][j] + H[0][i][1][j], 2) / (2. * variance));
					}
					else
					{
						app[m][0] *= exp(-pow(postRx[i][1][j] - estimate[0][i * SLIDING][j][1], 2) / (2. * variance));
						app[m][1] *= exp(-pow(postRx[i][1][j] + estimate[0][i * SLIDING][j][1], 2) / (2. * variance));
					}
				}
				else if (NUM_USER == 2)
				{
					if (CE_SCHEME == 1)
					{
						app[m][0] *= exp(-pow(postRx[i][1][j] - H[0][i][1][j] - H[1][i][1][j], 2) / (2. * variance));
						app[m][1] *= exp(-pow(postRx[i][1][j] - H[0][i][1][j] + H[1][i][1][j], 2) / (2. * variance));
						app[m][2] *= exp(-pow(postRx[i][1][j] + H[0][i][1][j] - H[1][i][1][j], 2) / (2. * variance));
						app[m][3] *= exp(-pow(postRx[i][1][j] + H[0][i][1][j] + H[1][i][1][j], 2) / (2. * variance));
					}
					else
					{
						app[m][0] *= exp(-pow(postRx[i][1][j] - estimate[0][i * SLIDING][j][1] - estimate[1][i * SLIDING][j][1], 2) / (2. * variance));
						app[m][1] *= exp(-pow(postRx[i][1][j] - estimate[0][i * SLIDING][j][1] + estimate[1][i * SLIDING][j][1], 2) / (2. * variance));
						app[m][2] *= exp(-pow(postRx[i][1][j] + estimate[0][i * SLIDING][j][1] - estimate[1][i * SLIDING][j][1], 2) / (2. * variance));
						app[m][3] *= exp(-pow(postRx[i][1][j] + estimate[0][i * SLIDING][j][1] + estimate[1][i * SLIDING][j][1], 2) / (2. * variance));
					}
				}
				else if (NUM_USER == 3)
				{
					if (CE_SCHEME == 1)
					{
						app[m][0] *= exp(-pow(postRx[i][1][j] - H[0][i][1][j] - H[1][i][1][j] - H[2][i][1][j], 2) / (2. * variance));
						app[m][1] *= exp(-pow(postRx[i][1][j] - H[0][i][1][j] - H[1][i][1][j] + H[2][i][1][j], 2) / (2. * variance));
						app[m][2] *= exp(-pow(postRx[i][1][j] - H[0][i][1][j] + H[1][i][1][j] - H[2][i][1][j], 2) / (2. * variance));
						app[m][3] *= exp(-pow(postRx[i][1][j] - H[0][i][1][j] + H[1][i][1][j] + H[2][i][1][j], 2) / (2. * variance));
						app[m][4] *= exp(-pow(postRx[i][1][j] + H[0][i][1][j] - H[1][i][1][j] - H[2][i][1][j], 2) / (2. * variance));
						app[m][5] *= exp(-pow(postRx[i][1][j] + H[0][i][1][j] - H[1][i][1][j] + H[2][i][1][j], 2) / (2. * variance));
						app[m][6] *= exp(-pow(postRx[i][1][j] + H[0][i][1][j] + H[1][i][1][j] - H[2][i][1][j], 2) / (2. * variance));
						app[m][7] *= exp(-pow(postRx[i][1][j] + H[0][i][1][j] + H[1][i][1][j] + H[2][i][1][j], 2) / (2. * variance));
					}
					else
					{
						app[m][0] *= exp(-pow(postRx[i][1][j] - estimate[0][i * SLIDING][j][1] - estimate[1][i * SLIDING][j][1] - estimate[2][i * SLIDING][j][1], 2) / (2. * variance));
						app[m][1] *= exp(-pow(postRx[i][1][j] - estimate[0][i * SLIDING][j][1] - estimate[1][i * SLIDING][j][1] + estimate[2][i * SLIDING][j][1], 2) / (2. * variance));
						app[m][2] *= exp(-pow(postRx[i][1][j] - estimate[0][i * SLIDING][j][1] + estimate[1][i * SLIDING][j][1] - estimate[2][i * SLIDING][j][1], 2) / (2. * variance));
						app[m][3] *= exp(-pow(postRx[i][1][j] - estimate[0][i * SLIDING][j][1] + estimate[1][i * SLIDING][j][1] + estimate[2][i * SLIDING][j][1], 2) / (2. * variance));
						app[m][4] *= exp(-pow(postRx[i][1][j] + estimate[0][i * SLIDING][j][1] - estimate[1][i * SLIDING][j][1] - estimate[2][i * SLIDING][j][1], 2) / (2. * variance));
						app[m][5] *= exp(-pow(postRx[i][1][j] + estimate[0][i * SLIDING][j][1] - estimate[1][i * SLIDING][j][1] + estimate[2][i * SLIDING][j][1], 2) / (2. * variance));
						app[m][6] *= exp(-pow(postRx[i][1][j] + estimate[0][i * SLIDING][j][1] + estimate[1][i * SLIDING][j][1] - estimate[2][i * SLIDING][j][1], 2) / (2. * variance));
						app[m][7] *= exp(-pow(postRx[i][1][j] + estimate[0][i * SLIDING][j][1] + estimate[1][i * SLIDING][j][1] + estimate[2][i * SLIDING][j][1], 2) / (2. * variance));
					}
				}
				//---------- normalization ----------
				temp = 0;
				for (int k = 0; k < NUM_LEVEL; k++)
				{
					if (app[m][k] < NUMERIC_LIMIT) app[m][k] = NUMERIC_LIMIT;
					temp += app[m][k];
				}
				for (int k = 0; k < NUM_LEVEL; k++)
				{
					app[m][k] /= temp;
					if (app[m][k] < NUMERIC_LIMIT) app[m][k] = NUMERIC_LIMIT;
				}
				m++;
			}
		}
	}
	if (CH_CODING_TYPE && JCD)
	{
		(DIFF_ENC && JOINT_DEC) ? ldpc.JointGSPA(app, JOINT_IT, INNER_IT) : ldpc.GSPA(app, LDPC_IT);
	}
	//---------- a-posteriori LLRs ----------
	if (CH_CODING_TYPE && DIFF_ENC && !JCD)
	{
		if (NUM_USER == 1)
		{
			for (int i = 0; i < FFT_POINT; i++)
			{
				if (app[i][0] <= NUMERIC_LIMIT) refLlr[0][i] = -LLR_LIMIT;
				else if (app[i][1] <= NUMERIC_LIMIT) refLlr[0][i] = LLR_LIMIT;
				else refLlr[0][i] = log(app[i][0] / app[i][1]);
			}
		}
		else if (NUM_USER == 2)
		{
			for (int i = 0; i < FFT_POINT; i++)
			{
				//---------- user-A ----------
				if ((app[i][0] + app[i][2]) <= NUMERIC_LIMIT) refLlr[0][i] = -LLR_LIMIT;
				else if ((app[i][1] + app[i][3]) <= NUMERIC_LIMIT) refLlr[0][i] = LLR_LIMIT;
				else refLlr[0][i] = log((app[i][0] + app[i][2]) / (app[i][1] + app[i][3]));
				//---------- user-B ----------
				if ((app[i][0] + app[i][1]) <= NUMERIC_LIMIT) refLlr[1][i] = -LLR_LIMIT;
				else if ((app[i][2] + app[i][3]) <= NUMERIC_LIMIT) refLlr[1][i] = LLR_LIMIT;
				else refLlr[1][i] = log((app[i][0] + app[i][1]) / (app[i][2] + app[i][3]));
			}
		}
		else if (NUM_USER == 3)
		{
			for (int i = 0; i < FFT_POINT; i++)
			{
				//---------- user-A ----------
				if ((app[i][0] + app[i][1] + app[i][2] + app[i][3]) <= NUMERIC_LIMIT) refLlr[0][i] = -LLR_LIMIT;
				else if ((app[i][4] + app[i][5] + app[i][6] + app[i][7]) <= NUMERIC_LIMIT) refLlr[0][i] = LLR_LIMIT;
				else refLlr[0][i] = log((app[i][0] + app[i][1] + app[i][2] + app[i][3]) / (app[i][4] + app[i][5] + app[i][6] + app[i][7]));
				//---------- user-B ----------
				if ((app[i][0] + app[i][1] + app[i][4] + app[i][5]) <= NUMERIC_LIMIT) refLlr[1][i] = -LLR_LIMIT;
				else if ((app[i][2] + app[i][3] + app[i][6] + app[i][7]) <= NUMERIC_LIMIT) refLlr[1][i] = LLR_LIMIT;
				else refLlr[1][i] = log((app[i][0] + app[i][1] + app[i][4] + app[i][5]) / (app[i][2] + app[i][3] + app[i][6] + app[i][7]));
				//---------- user-C ----------
				if ((app[i][0] + app[i][2] + app[i][4] + app[i][6]) <= NUMERIC_LIMIT) refLlr[2][i] = -LLR_LIMIT;
				else if ((app[i][1] + app[i][3] + app[i][5] + app[i][7]) <= NUMERIC_LIMIT) refLlr[2][i] = LLR_LIMIT;
				else refLlr[2][i] = log((app[i][0] + app[i][2] + app[i][4] + app[i][6]) / (app[i][1] + app[i][3] + app[i][5] + app[i][7]));
			}
		}
	}
	if (NUM_USER == 1)
	{
		for (int i = 0; i < CODE_LEN; i++)
		{
			if (app[i][0] <= NUMERIC_LIMIT) appLlr[0][i - DIFF_ENC*FFT_POINT] = -LLR_LIMIT;
			else if (app[i][1] <= NUMERIC_LIMIT) appLlr[0][i - DIFF_ENC*FFT_POINT] = LLR_LIMIT;
			else appLlr[0][i - DIFF_ENC*FFT_POINT] = log(app[i][0] / app[i][1]);
		}
	}
	else if (NUM_USER == 2)
	{
		for (int i = 0; i < CODE_LEN; i++)
		{
			//---------- user-A ----------
			if ((app[i][0] + app[i][1]) <= NUMERIC_LIMIT) appLlr[1][i - (DIFF_ENC - JOINT_DEC * JCD) * FFT_POINT] = -LLR_LIMIT;
			else if ((app[i][2] + app[i][3]) <= NUMERIC_LIMIT) appLlr[1][i - (DIFF_ENC - JOINT_DEC * JCD) * FFT_POINT] = LLR_LIMIT;
			else appLlr[0][i - (DIFF_ENC - JOINT_DEC * JCD) * FFT_POINT] = log((app[i][0] + app[i][1]) / (app[i][2] + app[i][3]));
			//---------- user-B ----------
			if ((app[i][0] + app[i][2]) <= NUMERIC_LIMIT) appLlr[0][i - (DIFF_ENC - JOINT_DEC * JCD) * FFT_POINT] = -LLR_LIMIT;
			else if ((app[i][1] + app[i][3]) <= NUMERIC_LIMIT) appLlr[0][i - (DIFF_ENC - JOINT_DEC * JCD) * FFT_POINT] = LLR_LIMIT;
			else appLlr[1][i - (DIFF_ENC - JOINT_DEC * JCD) * FFT_POINT] = log((app[i][0] + app[i][2]) / (app[i][1] + app[i][3]));
		}
	}
	else if (NUM_USER == 3)
	{
		for (int i = 0; i < CODE_LEN; i++)
		{
			//---------- user-A ----------
			if ((app[i][0] + app[i][1] + app[i][2] + app[i][3]) <= NUMERIC_LIMIT) appLlr[0][i - (DIFF_ENC - JOINT_DEC*JCD)*FFT_POINT] = -LLR_LIMIT;
			else if ((app[i][4] + app[i][5] + app[i][6] + app[i][7]) <= NUMERIC_LIMIT) appLlr[0][i - (DIFF_ENC - JOINT_DEC*JCD)*FFT_POINT] = LLR_LIMIT;
			else appLlr[0][i - (DIFF_ENC - JOINT_DEC*JCD)*FFT_POINT] = log((app[i][0] + app[i][1] + app[i][2] + app[i][3]) / (app[i][4] + app[i][5] + app[i][6] + app[i][7]));
			//---------- user-B ----------
			if ((app[i][0] + app[i][1] + app[i][4] + app[i][5]) <= NUMERIC_LIMIT) appLlr[1][i - (DIFF_ENC - JOINT_DEC*JCD)*FFT_POINT] = -LLR_LIMIT;
			else if ((app[i][2] + app[i][3] + app[i][6] + app[i][7]) <= NUMERIC_LIMIT) appLlr[1][i - (DIFF_ENC - JOINT_DEC*JCD)*FFT_POINT] = LLR_LIMIT;
			else appLlr[1][i - (DIFF_ENC - JOINT_DEC*JCD)*FFT_POINT] = log((app[i][0] + app[i][1] + app[i][4] + app[i][5]) / (app[i][2] + app[i][3] + app[i][6] + app[i][7]));
			//---------- user-C ----------
			if ((app[i][0] + app[i][2] + app[i][4] + app[i][6]) <= NUMERIC_LIMIT) appLlr[2][i - (DIFF_ENC - JOINT_DEC*JCD)*FFT_POINT] = -LLR_LIMIT;
			else if ((app[i][1] + app[i][3] + app[i][5] + app[i][7]) <= NUMERIC_LIMIT) appLlr[2][i - (DIFF_ENC - JOINT_DEC*JCD)*FFT_POINT] = LLR_LIMIT;
			else appLlr[2][i - (DIFF_ENC - JOINT_DEC*JCD)*FFT_POINT] = log((app[i][0] + app[i][2] + app[i][4] + app[i][6]) / (app[i][1] + app[i][3] + app[i][5] + app[i][7]));
		}
	}
}