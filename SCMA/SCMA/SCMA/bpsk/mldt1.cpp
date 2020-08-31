#define _USE_MATH_DEFINES
#include <iostream>
#include "parameters.h"
using namespace std;

void MLDT(double variance, double ***chCoef, double **rx, double **app, vector<vector<double>> &appLlr, double **estimate)
{
	int count = 0;
	for (int i = 0; i < CODE_LEN; i++)
	{
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
					if (CE_SCHEME == 1)
						estimate_sum += pow(-1, (reg % 2)) * chCoef[j][0][i / FADING_SIZE];
					else
						estimate_sum += pow(-1, (reg % 2)) * estimate[j][0];
					
					reg /= 2;
				}
				app[i][num_level] = exp(-pow(rx[i][0] - estimate_sum, 2) / (2. * (variance + count * MODIFIER)));
			}
				
			if (MODIFIED_LLR)
			{
				for (int j = 0; j < NUM_LEVEL; j++)
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
						estimate_sum += pow(-1, (reg % 2)) * chCoef[j][1][i / FADING_SIZE];
					
					else
						estimate_sum += pow(-1, (reg % 2)) * estimate[j][1];
					
					reg /= 2;
				}
				app[i][num_level] *= exp(-pow(rx[i][1] - estimate_sum, 2) / (2. * (variance + count * MODIFIER)));
			}
		
			if (MODIFIED_LLR)
			{
				for (int j = 0; j < NUM_LEVEL; j++)
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
	for (int i = 0; i < CODE_LEN; i++)
	{
		for (int nuser = 0; nuser < NUM_USER; nuser++)
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
			if (app_sum[0] <= NUMERIC_LIMIT) appLlr[nuser][i] = -LLR_LIMIT;
			else if (app_sum[1] <= NUMERIC_LIMIT) appLlr[nuser][i] = LLR_LIMIT;
			else appLlr[nuser][i] = log(app_sum[0] / app_sum[1]);
		}	
	}
}