#define _USE_MATH_DEFINES
#include <iostream>
#include "parameters.h"
#include <cmath>
using namespace std;

void MLDT(double variance, double **chCoef, double **rx, double **app, double **appLlr, double **estimate, int *known_drift, vector<vector<double>> map_16QAM)
{
	int count = 0;
	if (COLLISION)
	{
		for (int i = 0; i < BLOCK_LEN / MOD_LEVEL; i++)
		{//TDDOO
			while (true)
			{
				bool flag = true;
				//---------- APPs ----------
				for (int num_level = 0; num_level < NUM_LEVEL; num_level++)
				{
					int reg = num_level;
					double estimate_sum = 0;
					for (int j = NUM_USER - 1; j >= 0; j--)
					{
						if (CE_SCHEME == 1){
                            estimate_sum = map_16QAM[num_level][0] * chCoef[j][0] - map_16QAM[num_level][1] * chCoef[j][1];
						}
						else
							estimate_sum = map_16QAM[num_level][0] * estimate[j][0] - map_16QAM[num_level][1] * estimate[j][1];

						reg /= 2;
					}
					app[i][num_level] = exp(-pow(rx[i][0] - estimate_sum, 2) / (2. * (variance + count * MODIFIER)));
				}


				//---------- normalization ----------
				double temp = 0;
				for (int j = 0; j < NUM_LEVEL; j++)
				{
					if (app[i][j] < NUMERIC_LIMIT) app[i][j] = NUMERIC_LIMIT;
					temp += app[i][j];
				}
				for (int j = 0; j < NUM_LEVEL; j++)
				{
					app[i][j] /= temp;
					if (app[i][j] < NUMERIC_LIMIT) app[i][j] = NUMERIC_LIMIT;
				}
				//---------- APPs ----------
				for (int num_level = 0; num_level < NUM_LEVEL; num_level++)
				{
					int reg = num_level;
					double estimate_sum = 0;
					for (int j = NUM_USER - 1; j >= 0; j--)
					{
						if (CE_SCHEME == 1)
							estimate_sum = map_16QAM[num_level][0] * chCoef[j][1] + map_16QAM[num_level][1] * chCoef[j][0];

						else
							estimate_sum = map_16QAM[num_level][0] * estimate[j][1] + map_16QAM[num_level][1] * estimate[j][0];

						reg /= 2;
					}
					app[i][num_level] *= exp(-pow(rx[i][1] - estimate_sum, 2) / (2. * (variance + count * MODIFIER)));
				}

				//---------- normalization ----------
				temp = 0;
				for (int j = 0; j < NUM_LEVEL; j++)
				{
					if (app[i][j] < NUMERIC_LIMIT) app[i][j] = NUMERIC_LIMIT;
					temp += app[i][j];
				}
				for (int j = 0; j < NUM_LEVEL; j++)
				{
					app[i][j] /= temp;
					if (app[i][j] < NUMERIC_LIMIT) app[i][j] = NUMERIC_LIMIT;
				}


				if (flag) break;
			}
		}

		//---------- a-posteriori LLRs ----------
		for (int i = 0; i < BLOCK_LEN / MOD_LEVEL; i++)
		{
			for (int nuser = 0; nuser < NUM_USER*MOD_LEVEL; nuser++)
			{




				int range = NUM_LEVEL / pow(2, nuser + 1);
				double app_sum[2] = { 0 };
				for (int j = 0; j < NUM_LEVEL; j++)
				{
					if (j % (2 * range) < range)
						app_sum[0] += app[i][j];
					else
						app_sum[1] += app[i][j];
				}
				if (app_sum[0] <= NUMERIC_LIMIT) appLlr[0][4*i+nuser] = -LLR_LIMIT;
				else if (app_sum[1] <= NUMERIC_LIMIT) appLlr[0][4*i+nuser] = LLR_LIMIT;
				else appLlr[0][4*i+nuser] = log(app_sum[0] / app_sum[1]);
			}
		}
	}
	/*
	else
	{
		int R = 0, CORRENT_USER=0, CORRENT_LEVEL=0, step=0;
		// R is the bonder of each range
		// Step is parameter for using channel ceofficient


		for (int i = 0; i < known_drift[2 * NUM_USER - 1]; i++)
		{
			if (i >= known_drift[R])
			{
				R >= NUM_USER ? step++ : step ;
				R >= NUM_USER ? CORRENT_USER-- : CORRENT_USER++;
				CORRENT_LEVEL = pow(2, CORRENT_USER);
				R++;
			}
			
			while (true)
			{
				bool flag = true;
				//---------- APPs ----------
				for (int num_level = 0; num_level < CORRENT_LEVEL; num_level++)
				{
					int reg = num_level;
					double estimate_sum = 0;
					for (int j = CORRENT_USER - 1; j >= 0; j--)
					{
						if (CE_SCHEME == 1)
							estimate_sum += pow(-1, (reg % 2)) * chCoef[step + j][0];
						else
							estimate_sum += pow(-1, (reg % 2)) * estimate[step + j][0];

						reg /= 2;
					}
					app[i][num_level] = exp(-pow(rx[i][0] - estimate_sum, 2) / (2. * (variance + count * MODIFIER)));
				}

				if (MODIFIED_LLR)
				{
					for (int j = 0; j < CORRENT_LEVEL; j++)
					{
						if (app[i][j] == 0)
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
				for (int j = 0; j < CORRENT_LEVEL; j++)
				{
					if (app[i][j] < NUMERIC_LIMIT) app[i][j] = NUMERIC_LIMIT;
					temp += app[i][j];
				}
				for (int j = 0; j < CORRENT_LEVEL; j++)
				{
					app[i][j] /= temp;
					if (app[i][j] < NUMERIC_LIMIT) app[i][j] = NUMERIC_LIMIT;
				}
				//---------- APPs ----------
				for (int num_level = 0; num_level < CORRENT_LEVEL; num_level++)
				{
					int reg = num_level;
					double estimate_sum = 0;
					for (int j = CORRENT_USER - 1; j >= 0; j--)
					{
						if (CE_SCHEME == 1)
							estimate_sum += pow(-1, (reg % 2)) * chCoef[step + j][1];

						else
							estimate_sum += pow(-1, (reg % 2)) * estimate[step + j][1];

						reg /= 2;
					}
					app[i][num_level] *= exp(-pow(rx[i][1] - estimate_sum, 2) / (2. * (variance + count * MODIFIER)));
				}

				if (MODIFIED_LLR)
				{
					for (int j = 0; j < CORRENT_LEVEL; j++)
					{
						if (app[i][j] == 0)
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
				for (int j = 0; j < CORRENT_LEVEL; j++)
				{
					if (app[i][j] < NUMERIC_LIMIT) app[i][j] = NUMERIC_LIMIT;
					temp += app[i][j];
				}
				for (int j = 0; j < CORRENT_LEVEL; j++)
				{
					app[i][j] /= temp;
					if (app[i][j] < NUMERIC_LIMIT) app[i][j] = NUMERIC_LIMIT;
				}
				if (flag) break;
			}
			//---------- a-posteriori LLRs ----------

			for (int nuser = 0; nuser < CORRENT_USER; nuser++)
			{
				int range = CORRENT_LEVEL / pow(2, nuser + 1);
				double app_sum[2] = { 0 };
				for (int j = 0; j < CORRENT_LEVEL; j++)
				{
					if (j % (2 * range) < range)
						app_sum[0] += app[i][j];
					else
						app_sum[1] += app[i][j];
				}
				if (app_sum[0] <= NUMERIC_LIMIT) appLlr[nuser + step][i] = -LLR_LIMIT;
				else if (app_sum[1] <= NUMERIC_LIMIT) appLlr[nuser + step][i] = LLR_LIMIT;
				else appLlr[nuser + step][i] = log(app_sum[0] / app_sum[1]);
			}
		}
	}
	*/
}