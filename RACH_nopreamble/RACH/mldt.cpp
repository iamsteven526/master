#define _USE_MATH_DEFINES
#include <iostream>
#include "parameters.h"
#include <cstring>
#include <cmath>
#include <random>
using namespace std;

namespace
{
	random_device seed;
	mt19937 generator(seed());
	normal_distribution<double> normal(0, 1);
}

void MLDT(LDPC &ldpc, double stdDev, double *****H, double *****postRx, double ****app, double ***appLlr, double ***refLlr, double *****estimate, int* packet_num, int*** Cluster_num, int**** Cluster_gain, double snrdB)
{

	/*for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		for (int t = 0; t < packet_num[nuser]; t++)
		{
			for (int i = 0, m = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
			{
				for (int j = 0; j < FFT_POINT; j++)
				{
					cout << estimate[0][0][i * SLIDING][j][0]<<" "<<H[0][0][i][0][j]<< endl;

					if (abs(estimate[0][0][i * SLIDING][j][0]) < 0.0001 || abs(estimate[0][0][i * SLIDING][j][1]) < 0.0001)
					{
						cout << estimate[0][0][i * SLIDING][j][0] << " " << estimate[0][0][i * SLIDING][j][1] << endl;
						system("pause");
					}
				}
				
			}
		}
	}
	system("pause");*/
	//---- Sperate MLDT ---- 
	double adapted_variance = pow(stdDev + (0.1 * sqrt(snrdB)), 2);
	double test = 0,testtime = 0;
	double alpha = 0;
	int NUM_user;
	int NUM_level;
	//double adapted_variance = pow(stdDev, 2);
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		for (int t = 0; t < packet_num[nuser]; t++)
		{
			//0925modify
            for (int i = 0, m = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
			{
				//---- u_th user, t_th packet, i-th subframe ---- 
				NUM_user = Cluster_num[t][nuser][i];
				NUM_level = pow(2, Cluster_num[t][nuser][i]);
				
				//---- APPs_real_part ----
				for (int j = 0; j < FFT_POINT; j++)
				{
					for (int k = 0; k < NUM_level; k++)
					{
						int reg = k;
						double estimate_sum = 0;
						for (int g = NUM_user - 1; g >= 0; g--)//---- g is the variable to selete special channel gain "Cluster_gain[t][nuser][i][g]"
						{
							estimate_sum += pow(-1, (reg % 2)) * H[0][Cluster_gain[t][nuser][i][g]][i][0][j];

							reg /= 2;
						}
						test += pow(postRx[t][nuser][i][0][j] - estimate_sum, 2);
						//cout << pow(postRx[t][nuser][i][0][j] - estimate_sum, 2) << endl;
						testtime += 1;
					}
					for (int k = 0; k < NUM_level; k++)
					{
						int reg = k;
						double estimate_sum = 0;
						for (int g = NUM_user - 1; g >= 0; g--)
						{
							estimate_sum += pow(-1, (reg % 2)) * H[0][Cluster_gain[t][nuser][i][g]][i][1][j];

							reg /= 2;
						}
						test += pow(postRx[t][nuser][i][0][j] - estimate_sum, 2);
					}
					m++;					
				}
			}

			//0925endmodify
			//cout << testtime << endl;
			/*CSI
			if(testtime > 32500) alpha = 1;//user >= 6
			else if(testtime > 16400) alpha = 0.9;//user = 5
			else if(testtime > 8300) alpha = 0.8;//user = 4
			else if(testtime > 4200) alpha = 0.6;//user = 3
			else if(testtime < 2200) alpha = 0.4;//user = 1
			else alpha = 0.4;//user = 2
			*/
			if(testtime > 32500) alpha = 1;//user >= 6
			else if(testtime > 16400) alpha = 0.9;//user = 5
			else if(testtime > 8300) alpha = 0.8;//user = 4
			else if(testtime > 4200) alpha = 0.6;//user = 3
			else if(testtime < 2200) alpha = 0.4;//user = 1
			else alpha = 0.4;//user = 2		
            adapted_variance = sqrt(test) / (alpha*testtime);


			for (int i = 0, m = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
			{
				//---- u_th user, t_th packet, i-th subframe ---- 
				NUM_user = Cluster_num[t][nuser][i];
				NUM_level = pow(2, Cluster_num[t][nuser][i]);
				
				//---- APPs_real_part ----
				for (int j = 0; j < FFT_POINT; j++)
				{
					for (int k = 0; k < NUM_level; k++)
					{
						int reg = k;
						double estimate_sum = 0;
						for (int g = NUM_user - 1; g >= 0; g--)//---- g is the variable to selete special channel gain "Cluster_gain[t][nuser][i][g]"
						{
							if (CE_SCHEME == 1)
								estimate_sum += pow(-1, (reg % 2)) * H[0][Cluster_gain[t][nuser][i][g]][i][0][j];
							else
								//estimate_sum += pow(-1, (reg % 2)) * (estimate[0][Cluster_gain[t][nuser][i][g]][i * SLIDING][j][0] + 0.2*(pow(NUM_user,8)/600) * sqrt(0.01275 / FFT_SEGMENT) * normal(generator));
								estimate_sum += pow(-1, (reg % 2)) * (H[0][Cluster_gain[t][nuser][i][g]][i][0][j] + 0.2*(pow(NUM_user,8)/600) * sqrt(0.01275 / FFT_SEGMENT) * normal(generator));
							reg /= 2;
						}
						app[t][nuser][m][k] = exp(-pow(postRx[t][nuser][i][0][j] - estimate_sum, 2) / (2. * adapted_variance));
						//cout << app[t][nuser][m][k] << " ";
					}
					//system("pause");

					
					double temp = 0;
					for (int k = 0; k < NUM_level; k++)
					{
						if (app[t][nuser][m][k] < NUMERIC_LIMIT) app[t][nuser][m][k] = NUMERIC_LIMIT;
						temp += app[t][nuser][m][k];
					}
					for (int k = 0; k < NUM_level; k++)
					{
						app[t][nuser][m][k] /= temp;
						if (app[t][nuser][m][k] < NUMERIC_LIMIT) app[t][nuser][m][k] = NUMERIC_LIMIT;
					}

					for (int k = 0; k < NUM_level; k++)
					{
						int reg = k;
						double estimate_sum = 0;
						for (int g = NUM_user - 1; g >= 0; g--)
						{
							if (CE_SCHEME == 1)
								estimate_sum += pow(-1, (reg % 2)) * H[0][Cluster_gain[t][nuser][i][g]][i][1][j];
							else
								//estimate_sum+= pow(-1, (reg % 2)) * (estimate[0][Cluster_gain[t][nuser][i][g]][i * SLIDING][j][1] +  0.2*(pow(NUM_user,8)/600) * sqrt(0.01275 / FFT_SEGMENT) * normal(generator));
								estimate_sum+= pow(-1, (reg % 2)) * (H[0][Cluster_gain[t][nuser][i][g]][i][1][j] +  0.2*(pow(NUM_user,8)/600) * sqrt(0.01275 / FFT_SEGMENT) * normal(generator));
							reg /= 2;
						}
						app[t][nuser][m][k] *= exp(-pow(postRx[t][nuser][i][1][j] - estimate_sum, 2) / (2. * adapted_variance));

					//	cout << app[t][nuser][m][k] << " ";
					}
				
						

						

					//---------- normalization ----------
					temp = 0;
					for (int k = 0; k < NUM_level; k++)
					{
						if (app[t][nuser][m][k] < NUMERIC_LIMIT) app[t][nuser][m][k] = NUMERIC_LIMIT;
						temp += app[t][nuser][m][k];
					}
					for (int k = 0; k < NUM_level; k++)
					{
						app[t][nuser][m][k] /= temp;
						if (app[t][nuser][m][k] < NUMERIC_LIMIT) app[t][nuser][m][k] = NUMERIC_LIMIT;
						//cout << app[t][nuser][m][k] << " ";
					}
					//system("pause");
					m++;

					
					
				}
			}

			//---------- a-posteriori LLRs ----------

			//---- refLlr ----
			//---- diff_enc_llr ----
			if (DIFF_ENC && !JCD)
			{
				int NUM_user = Cluster_num[t][nuser][0];
				int NUM_level = pow(2, Cluster_num[t][nuser][0]);
				int ID_user;
				for (int x = 0; x < NUM_user; x++)
					if (nuser == Cluster_gain[t][nuser][0][x])
						ID_user = x;
				int range = NUM_level / pow(2, ID_user + 1);
				
				for (int i = 0; i < FFT_POINT; i++)
				{
					double app_sum[2] = { 0 };
					for (int j = 0; j < NUM_level; j++)
					{
						if (j % (2 * range) < range)
							app_sum[0] += app[t][nuser][i][j];
						else
							app_sum[1] += app[t][nuser][i][j];
					}
					
					if (app_sum[0] <= NUMERIC_LIMIT) refLlr[t][nuser][i] = -LLR_LIMIT;
					else if (app_sum[1] <= NUMERIC_LIMIT) refLlr[t][nuser][i] = LLR_LIMIT;
					else refLlr[t][nuser][i] = log(app_sum[0] / app_sum[1]);
				}	
			}

			//---- appLlr ----
			for (int i = DIFF_ENC * FFT_POINT; i < CODE_LEN + DIFF_ENC * FFT_POINT; i++)
			{
				int NUM_user = Cluster_num[t][nuser][i / FFT_POINT];
				int NUM_level = pow(2, Cluster_num[t][nuser][i / FFT_POINT]);
				int ID_user; 
				for (int x = 0; x < NUM_user; x++) 
					if (nuser == Cluster_gain[t][nuser][i / FFT_POINT][x]) 
						  ID_user = x;  
				int range = NUM_level / pow(2, ID_user + 1);
				double app_sum[2] = { 0 };

				for (int j = 0; j < NUM_level; j++)
				{
					if (j % (2 * range) < range)
						app_sum[0] += app[t][nuser][i][j];
					else
						app_sum[1] += app[t][nuser][i][j];
				}

				if (app_sum[0] <= NUMERIC_LIMIT) appLlr[t][nuser][i - DIFF_ENC * FFT_POINT] = -LLR_LIMIT;
				else if (app_sum[1] <= NUMERIC_LIMIT) appLlr[t][nuser][i - DIFF_ENC * FFT_POINT] = LLR_LIMIT;
				else appLlr[t][nuser][i - DIFF_ENC * FFT_POINT] = log(app_sum[0] / app_sum[1]);

				//cout << app_sum[0] << " "<< app_sum[1];
				//system("pause");
				//appLlr[t][nuser][i - DIFF_ENC * FFT_POINT] = appLlr[t][nuser][i - DIFF_ENC * FFT_POINT] > LLR_LIMIT ? LLR_LIMIT : appLlr[t][nuser][i - DIFF_ENC * FFT_POINT];
				//appLlr[t][nuser][i - DIFF_ENC * FFT_POINT] = appLlr[t][nuser][i - DIFF_ENC * FFT_POINT] < -LLR_LIMIT ? -LLR_LIMIT : appLlr[t][nuser][i - DIFF_ENC * FFT_POINT];
			}
			//system("pause");
			test = 0;
			testtime = 0;
		}
	}

	/*for (int i = 0; i < CODE_LEN; i++)
		cout << appLlr[0][0][i] << " ";

	system("pause");*/
	//cout << "est : " << NUM_user << endl;
}
