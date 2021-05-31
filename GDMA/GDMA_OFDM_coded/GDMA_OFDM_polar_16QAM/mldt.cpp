#define _USE_MATH_DEFINES
#include <iostream>
#include "parameters.h"
#include <random>
using namespace std;

namespace
{
	random_device seed;
	mt19937 generator(seed());
	normal_distribution<double> normal(0, 1);
}

void MLDT(LDPC &ldpc, double variance, double ****H, double ***postRx, double **app, double **appLlr, double **refLlr, double ****estimate, vector<vector<double>> map_16QAM)
{
	
	/*for (int j = 0; j < FFT_POINT; j++)
	{
		int temp = rand() % 2;
		for (int i = 0; i < FFT_SEGMENT; i++)
		{
			H[0][i][0][j] *= (2 * temp - 1);
			H[0][i][1][j] *= (2 * temp - 1);

			//cout << H[0][i][1][j] << " ";
		}
		//cout << endl;
	}*/
	//system("pause");
	if (CE_SCHEME == 0){
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			for (int j = 0; j < FFT_POINT; j++)
			{
				int temp = 1;//rand() % 2;
				for (int i = 0; i < FFT_SEGMENT; i++)
				{
					estimate[nuser][i][j][0] = H[nuser][i][0][j] * (2 * temp - 1) + sqrt(1.5*variance / (4*FFT_SEGMENT)) * normal(generator);
					estimate[nuser][i][j][1] = H[nuser][i][1][j] * (2 * temp - 1) + sqrt(1.5*variance / (4*FFT_SEGMENT)) * normal(generator);
				}
			}
		}
	}
	for (int i = 0, m = 0; i < FFT_SEGMENT; i++)
	{
		for (int j = 0; j < FFT_POINT; j++)
		{
			//---------- APPs ----------
			if (NUM_USER == 1)
			{
				if (CE_SCHEME == 1)
				{
					for(int k = 0; k < NUM_LEVEL; k++){
						app[m][k] = exp(-pow(postRx[i][0][j] - (map_16QAM[k][0]*H[0][i][0][j] - map_16QAM[k][1]*H[0][i][1][j]), 2) / (2. * variance));
					}
					//app[m][0] = exp(-pow(postRx[i][0][j] - H[0][i][0][j], 2) / (2. * variance));
					//app[m][1] = exp(-pow(postRx[i][0][j] + H[0][i][0][j], 2) / (2. * variance));
				}
				else
				{
					for(int k = 0; k < NUM_LEVEL; k++){
						app[m][k] = exp(-pow(postRx[i][0][j] - (map_16QAM[k][0]*estimate[0][i][j][0] - map_16QAM[k][1]*estimate[0][i][j][1]), 2) / (2. * variance));
					}
					//app[m][0] = exp(-pow(postRx[i][0][j] - estimate[0][i * SLIDING][j][0], 2) / (2. * variance));
					//app[m][1] = exp(-pow(postRx[i][0][j] + estimate[0][i * SLIDING][j][0], 2) / (2. * variance));
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
				//system("pause");
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
					for(int k = 0; k < NUM_LEVEL; k++){
						app[m][k] *= exp(-pow(postRx[i][1][j] - (map_16QAM[k][0]*H[0][i][1][j] + map_16QAM[k][1]*H[0][i][0][j]), 2) / (2. * variance));
					}
					//app[m][0] *= exp(-pow(postRx[i][1][j] - H[0][i][1][j], 2) / (2. * variance));
					//app[m][1] *= exp(-pow(postRx[i][1][j] + H[0][i][1][j], 2) / (2. * variance));
				}
				else
				{
					for(int k = 0; k < NUM_LEVEL; k++){
						app[m][k] *= exp(-pow(postRx[i][1][j] - (map_16QAM[k][0]*estimate[0][i][j][1] + map_16QAM[k][1]*estimate[0][i][j][0]), 2) / (2. * variance));
					}
					//app[m][0] *= exp(-pow(postRx[i][1][j] - estimate[0][i * SLIDING][j][1], 2) / (2. * variance));
					//app[m][1] *= exp(-pow(postRx[i][1][j] + estimate[0][i * SLIDING][j][1], 2) / (2. * variance));
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
	
    //LLR part
	if (NUM_USER == 1)
	{
		for (int i = 0; i < CODE_LEN / MOD_LEVEL; i++)
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



			//if (app[i][0] <= NUMERIC_LIMIT) appLlr[0][i] = -LLR_LIMIT;
			//else if (app[i][1] <= NUMERIC_LIMIT) appLlr[0][i] = LLR_LIMIT;
			//else appLlr[0][i] = log(app[i][0] / app[i][1]);
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

	/*for (int i = 0; i < FFT_POINT; i++)
	{
		for (int j = 0; j < FFT_SEGMENT; j++)
		{
			cout << HARD(appLlr[0][j*FFT_POINT+i]) << " ";
		}
		cout << endl;
	}*/

	/*for (int i = 0; i < CODE_LEN; i++)
	{
		cout << appLlr[0][i] << " ";
	}
	system("pause");
	*/
}