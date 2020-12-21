#include <iostream>
#include "parameters.h"
#include <cstring>
#include <math.h>
#include "scma_log_mpa.h"
#include <omp.h>
using namespace std;

void AlamoutiEncoder(int **data, double ***tx)
{
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		//---------- messages generation ----------
		data[nuser][0] = 1;//rand() % 2;
		data[nuser][1] = 0;//rand() % 2;
		//---------- Alamouti encoding ----------
		tx[nuser][0][0] = 1 - 2 * data[nuser][0];
		tx[nuser][1][0] = 1 - 2 * data[nuser][1];
		tx[nuser][0][1] = -tx[nuser][1][0];
		tx[nuser][1][1] = +tx[nuser][0][0];
	}
}

void SignalCombiner(double ***chCoef, double **rx, double ***postRx)
{
	for (int i = 0; i < NUM_USER; i++)
	{
		//double nFactor = sqrt(pow(chCoef[i][0][0], 2) + pow(chCoef[i][0][1], 2) + pow(chCoef[i][1][0], 2) + pow(chCoef[i][1][1], 2));
		double nFactor = 1;
		postRx[i][0][0] = (chCoef[i][0][0] * rx[0][0] + chCoef[i][0][1] * rx[0][1] + chCoef[i][1][0] * rx[1][0] + chCoef[i][1][1] * rx[1][1]) / nFactor;
		postRx[i][0][1] = (chCoef[i][0][0] * rx[0][1] - chCoef[i][0][1] * rx[0][0] + chCoef[i][1][1] * rx[1][0] - chCoef[i][1][0] * rx[1][1]) / nFactor;
		postRx[i][1][0] = (chCoef[i][1][0] * rx[0][0] + chCoef[i][1][1] * rx[0][1] - chCoef[i][0][0] * rx[1][0] - chCoef[i][0][1] * rx[1][1]) / nFactor;
		postRx[i][1][1] = (chCoef[i][1][0] * rx[0][1] - chCoef[i][1][1] * rx[0][0] - chCoef[i][0][1] * rx[1][0] + chCoef[i][0][0] * rx[1][1]) / nFactor;
	}
}

void ComplexConjugate(double *x, double *y)
{
	y[0] = x[0]; y[1] = -x[1];
}

void ComplexMultiplication(double *x1, double *x2, double *y)
{
	double temp[2];
	temp[0] = x1[0] * x2[0] - x1[1] * x2[1];
	temp[1] = x1[1] * x2[0] + x1[0] * x2[1];
	y[0] = temp[0]; y[1] = temp[1];
}

void SuperLevelSpecification(double ***chCoef, double ****supLevel)
{
	if (NUM_USER == 1)
	{
		supLevel[0][0][0][0] = supLevel[0][1][0][0] = pow(chCoef[0][0][0], 2) + pow(chCoef[0][0][1], 2) + pow(chCoef[0][1][0], 2) + pow(chCoef[0][1][1], 2);
		supLevel[0][0][0][1] = supLevel[0][1][0][1] = 0;
	}
	else if (NUM_USER == 2)
	{
		double temp1[2], temp2[2];
		supLevel[0][0][0][0] = supLevel[0][1][0][0] = pow(chCoef[0][0][0], 2) + pow(chCoef[0][0][1], 2) + pow(chCoef[0][1][0], 2) + pow(chCoef[0][1][1], 2);
		supLevel[0][0][0][1] = supLevel[0][1][0][1] = 0;
		ComplexConjugate(chCoef[0][0], temp1);
		ComplexMultiplication(temp1, chCoef[1][0], temp1);
		ComplexConjugate(chCoef[1][1], temp2);
		ComplexMultiplication(temp2, chCoef[0][1], temp2);
		supLevel[0][0][1][0] = temp1[0] + temp2[0];
		supLevel[0][0][1][1] = temp1[1] + temp2[1];
		supLevel[1][1][2][0] = supLevel[0][0][1][0];
		supLevel[1][1][2][1] = supLevel[0][0][1][1];
		ComplexConjugate(supLevel[0][0][1], supLevel[0][1][2]);
		supLevel[1][0][1][0] = supLevel[0][1][2][0];
		supLevel[1][0][1][1] = supLevel[0][1][2][1];
		supLevel[1][0][0][0] = supLevel[1][1][0][0] = pow(chCoef[1][0][0], 2) + pow(chCoef[1][0][1], 2) + pow(chCoef[1][1][0], 2) + pow(chCoef[1][1][1], 2);
		supLevel[1][0][0][1] = supLevel[1][1][0][1] = 0;
		ComplexConjugate(chCoef[0][0], temp1);
		ComplexMultiplication(temp1, chCoef[1][1], temp1);
		ComplexConjugate(chCoef[1][0], temp2);
		ComplexMultiplication(temp2, chCoef[0][1], temp2);
		supLevel[0][0][2][0] = temp1[0] - temp2[0];
		supLevel[0][0][2][1] = temp1[1] - temp2[1];
		ComplexConjugate(supLevel[0][0][2], supLevel[1][1][1]);
		supLevel[1][0][2][0] = -supLevel[0][0][2][0];
		supLevel[1][0][2][1] = -supLevel[0][0][2][1];
		ComplexConjugate(supLevel[1][0][2], supLevel[0][1][1]);
	}
	else if (NUM_USER == 3)
	{
	    //superlevel[user][tx][user+tx-1][2]
		double temp1[2], temp2[2];	
		supLevel[0][0][0][0] = supLevel[0][1][0][0] = pow(chCoef[0][0][0], 2) + pow(chCoef[0][0][1], 2) + pow(chCoef[0][1][0], 2) + pow(chCoef[0][1][1], 2);
		supLevel[0][0][0][1] = supLevel[0][1][0][1] = 0;
		ComplexConjugate(chCoef[0][0], temp1);
		ComplexMultiplication(temp1, chCoef[1][0], temp1);
		ComplexConjugate(chCoef[1][1], temp2);
		ComplexMultiplication(temp2, chCoef[0][1], temp2);
		supLevel[0][0][1][0] = temp1[0] + temp2[0];
		supLevel[0][0][1][1] = temp1[1] + temp2[1];

		ComplexConjugate(chCoef[0][0], temp1);
		ComplexMultiplication(temp1, chCoef[2][0], temp1);
		ComplexConjugate(chCoef[2][1], temp2);
		ComplexMultiplication(temp2, chCoef[0][1], temp2);
		supLevel[0][0][3][0] = temp1[0] + temp2[0];
		supLevel[0][0][3][1] = temp1[1] + temp2[1];

		supLevel[1][1][2][0] = supLevel[0][0][1][0];
		supLevel[1][1][2][1] = supLevel[0][0][1][1];
		supLevel[2][1][2][0] = supLevel[0][0][3][0]; 
		supLevel[2][1][2][1] = supLevel[0][0][3][1]; 
		ComplexConjugate(supLevel[0][0][1], supLevel[0][1][2]);
		ComplexConjugate(supLevel[0][0][3], supLevel[0][1][4]);
		supLevel[1][0][1][0] = supLevel[0][1][2][0];
		supLevel[1][0][1][1] = supLevel[0][1][2][1];
		supLevel[2][0][1][0] = supLevel[0][1][4][0];
		supLevel[2][0][1][1] = supLevel[0][1][4][1];


		supLevel[1][0][0][0] = supLevel[1][1][0][0] = pow(chCoef[1][0][0], 2) + pow(chCoef[1][0][1], 2) + pow(chCoef[1][1][0], 2) + pow(chCoef[1][1][1], 2);
		supLevel[1][0][0][1] = supLevel[1][1][0][1] = 0;
		ComplexConjugate(chCoef[0][0], temp1);
		ComplexMultiplication(temp1, chCoef[1][1], temp1);
		ComplexConjugate(chCoef[1][0], temp2);
		ComplexMultiplication(temp2, chCoef[0][1], temp2);
		supLevel[0][0][2][0] = temp1[0] - temp2[0];
		supLevel[0][0][2][1] = temp1[1] - temp2[1];
		ComplexConjugate(supLevel[0][0][2], supLevel[1][1][1]);
		supLevel[1][0][2][0] = -supLevel[0][0][2][0];
		supLevel[1][0][2][1] = -supLevel[0][0][2][1];
		ComplexConjugate(supLevel[1][0][2], supLevel[0][1][1]);

		ComplexConjugate(chCoef[0][0], temp1);
		ComplexMultiplication(temp1, chCoef[2][1], temp1);
		ComplexConjugate(chCoef[2][0], temp2);
		ComplexMultiplication(temp2, chCoef[0][1], temp2);
		supLevel[0][0][4][0] = temp1[0] - temp2[0];
		supLevel[0][0][4][1] = temp1[1] - temp2[1];
		ComplexConjugate(supLevel[0][0][4], supLevel[2][1][1]);
		supLevel[2][0][2][0] = -supLevel[0][0][4][0];
		supLevel[2][0][2][1] = -supLevel[0][0][4][1];
		ComplexConjugate(supLevel[2][0][2], supLevel[0][1][3]);

		supLevel[2][0][0][0] = supLevel[2][1][0][0] = pow(chCoef[2][0][0], 2) + pow(chCoef[2][0][1], 2) + pow(chCoef[2][1][0], 2) + pow(chCoef[2][1][1], 2);
		supLevel[2][0][0][1] = supLevel[2][1][0][1] = 0;
		ComplexConjugate(chCoef[1][0], temp1);
		ComplexMultiplication(temp1, chCoef[2][0], temp1);
		ComplexConjugate(chCoef[2][1], temp2);
		ComplexMultiplication(temp2, chCoef[1][1], temp2);
		supLevel[1][0][3][0] = temp1[0] + temp2[0];
		supLevel[1][0][3][1] = temp1[1] + temp2[1];
		supLevel[2][1][4][0] = supLevel[1][0][3][0];
		supLevel[2][1][4][1] = supLevel[1][0][3][1];
		ComplexConjugate(supLevel[1][0][3], supLevel[1][1][4]);
		supLevel[2][0][3][0] = supLevel[1][1][4][0];
		supLevel[2][0][3][1] = supLevel[1][1][4][1];

		ComplexConjugate(chCoef[1][0], temp1);
		ComplexMultiplication(temp1, chCoef[2][1], temp1);
		ComplexConjugate(chCoef[2][0], temp2);
		ComplexMultiplication(temp2, chCoef[1][1], temp2);
		supLevel[0][1][4][0] = temp1[0] - temp2[0];
		supLevel[0][1][4][1] = temp1[1] - temp2[1];
		ComplexConjugate(supLevel[0][1][4], supLevel[2][1][3]);
		supLevel[2][0][4][0] = -supLevel[0][1][4][0];
		supLevel[2][0][4][1] = -supLevel[0][1][4][1];
		ComplexConjugate(supLevel[2][0][4], supLevel[1][1][3]);


	}
	else
	{
		printf("\nPARAMETER SETTING IS WRONG\n");
	}
	for (int i = 0; i < NUM_USER; i++)
	{
		//double nFactor = sqrt(pow(chCoef[i][0][0], 2) + pow(chCoef[i][0][1], 2) + pow(chCoef[i][1][0], 2) + pow(chCoef[i][1][1], 2));
		double nFactor = 1;
		for (int j = 0; j < NUM_TX; j++)
		{
			for (int k = 0; k < (NUM_USER*NUM_TX - 1); k++)
			{
				supLevel[i][j][k][0] /= nFactor;
				supLevel[i][j][k][1] /= nFactor;
			}
		}
	}
}

void Detector(int **data, double **appLlr, long double &error)
{
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		for (int i = 0; i < NUM_TX; i++)
		{			
			error += (data[nuser][i] != HARD(appLlr[nuser][i]));
		}
	}
}

void    CALC_F(int data_idx[SCMA_SOURCE][SCMA_SOURCE * NUM_USER / SCMA_USER_SOURCE], double ***appLlr)
{
	long double f[NUM_TX][K][M][M][M] = { 0 };
	long double m1p[M] = { 0 };
	long double m2p[M] = { 0 };
	long double m3p[M] = { 0 };
	for(int ntx = 0; ntx < NUM_TX; ntx++)
	{
		for (int k = 0; k < K; k++)
		{
			for (int m1 = 0; m1 < M; m1++)
			{
                m1p[0] = exp(appLlr[k][0][ntx])/(1.0 + exp(appLlr[k][0][ntx]));
				m1p[1] = 1.0/(1.0 + exp(appLlr[k][0][ntx]));
				//cout << m1p[0] << " " << m1p[1] << endl;
				for (int m2 = 0; m2 < M; m2++)
				{
	                m2p[0] = exp(appLlr[k][1][ntx])/(1.0 + exp(appLlr[k][1][ntx]));
				    m2p[1] = 1.0/(1.0 + exp(appLlr[k][1][ntx]));				
					for (int m3 = 0; m3 < M; m3++)
					{
		                m3p[0] = exp(appLlr[k][2][ntx])/(1.0 + exp(appLlr[k][2][ntx]));
				        m3p[1] = 1.0/(1.0 + exp(appLlr[k][2][ntx]));
						f[ntx][k][m1][m2][m3] = m1p[m1]*m2p[m2]*m3p[m3];
					}
				}
			}
		}
	}
	///////start tx0 bit
	double Ap = log(1.0 / M);
	double Igv[K][V][M] = { 0 };
	double Ivg[K][V][M] = { 0 };
	int Niter = 5;
	for (int k = 0; k < K; k++)
	{
		for (int v = 0; v < V; v++)
		{
			for (int m = 0; m < M; m++)
			{
				Ivg[k][v][m] = Ap;
			}
		}
	}
	int ind_df[SCMA_SOURCE][NUM_USER] = { {1,2,4}, {0,2,5}, {1,3,5}, {0,3,4} };
	int ind_dv[SCMA_SOURCE * NUM_USER / SCMA_USER_SOURCE][SCMA_USER_SOURCE] = { {1,3}, {0,2}, {0,1}, {2,3}, {0,3}, {1,2} };

	//iteration for tx0
	for (int iter = 0; iter < Niter; iter++)
	{
		// Igv update

		for (int k = 0; k < K; k++)
		{

			for (int m1 = 0; m1 < M; m1++)
			{
				double sIgv[M * M] = { 0 };
				for (int m2 = 0; m2 < M; m2++)
				{
					for (int m3 = 0; m3 < M; m3++)
					{
						sIgv[m2 * M + m3] = f[0][k][m1][m2][m3] + Ivg[k][ind_df[k][1]][m2] + Ivg[k][ind_df[k][2]][m3];
					}
				}
				Igv[k][ind_df[k][0]][m1] = log_sum_exp(sIgv, M * M);
			}

			for (int m2 = 0; m2 < M; m2++)
			{
				double sIgv[M * M] = { 0 };
				for (int m1 = 0; m1 < M; m1++)
				{
					for (int m3 = 0; m3 < M; m3++)
					{
						sIgv[m1 * M + m3] = f[0][k][m1][m2][m3] + Ivg[k][ind_df[k][0]][m1] + Ivg[k][ind_df[k][2]][m3];
					}
				}
				Igv[k][ind_df[k][1]][m2] = log_sum_exp(sIgv, M * M);
			}

			for (int m3 = 0; m3 < M; m3++)
			{
				double sIgv[M * M] = { 0 };
				for (int m1 = 0; m1 < M; m1++)
				{
					for (int m2 = 0; m2 < M; m2++)
					{
						sIgv[m1 * M + m2] = f[0][k][m1][m2][m3] + Ivg[k][ind_df[k][0]][m1] + Ivg[k][ind_df[k][1]][m2];
					}
				}
				Igv[k][ind_df[k][2]][m3] = log_sum_exp(sIgv, M * M);
			}
		}

		// Ivg update

		for (int v = 0; v < V; v++)
		{
			double sum0 = 0;
			double sum1 = 0;

			for (int i = 0; i < M; i++)
			{
				sum0 += exp(Igv[ind_dv[v][0]][v][i]);
				sum1 += exp(Igv[ind_dv[v][1]][v][i]);
			}

			sum0 = log(sum0);
			sum1 = log(sum1);

			for (int m = 0; m < M; m++)
			{
				Ivg[ind_dv[v][0]][v][m] = Igv[ind_dv[v][1]][v][m] - sum1;
				Ivg[ind_dv[v][1]][v][m] = Igv[ind_dv[v][0]][v][m] - sum0;
			}
		}
	}

	// Step 3: LLR calculation//TODO

	double Q[M][V] = { 0 };

	for (int v = 0; v < V; v++)
	{
		for (int m = 0; m < M; m++)
		{
			Q[m][v] = Ap + Igv[ind_dv[v][0]][v][m] + Igv[ind_dv[v][1]][v][m];
		}
	}

    for(int nresource = 0; nresource < K; nresource++)
	{	
		appLlr[nresource][0][0] = Q[0][ind_df[nresource][0]] - Q[1][ind_df[nresource][0]];
		appLlr[nresource][1][0] = Q[0][ind_df[nresource][1]] - Q[1][ind_df[nresource][1]];
		appLlr[nresource][2][0] = Q[0][ind_df[nresource][2]] - Q[1][ind_df[nresource][2]];
	}
	///////start tx1 bit
	Ap = log(1.0 / M);
	Igv[K][V][M] = { 0 };
	Ivg[K][V][M] = { 0 };
	Niter = 5;
	for (int k = 0; k < K; k++)
	{
		for (int v = 0; v < V; v++)
		{
			for (int m = 0; m < M; m++)
			{
				Ivg[k][v][m] = Ap;
			}
		}
	}
	//ind_df[SCMA_SOURCE][NUM_USER] = { {1,2,4}, {0,2,5}, {1,3,5}, {0,3,4} };
	//ind_dv[SCMA_SOURCE * NUM_USER / SCMA_USER_SOURCE][SCMA_USER_SOURCE] = { {1,3}, {0,2}, {0,1}, {2,3}, {0,3}, {1,2} };

	//iteration for tx0
	for (int iter = 0; iter < Niter; iter++)
	{
		// Igv update

		for (int k = 0; k < K; k++)
		{

			for (int m1 = 0; m1 < M; m1++)
			{
				double sIgv[M * M] = { 0 };
				for (int m2 = 0; m2 < M; m2++)
				{
					for (int m3 = 0; m3 < M; m3++)
					{
						sIgv[m2 * M + m3] = f[1][k][m1][m2][m3] + Ivg[k][ind_df[k][1]][m2] + Ivg[k][ind_df[k][2]][m3];
					}
				}
				Igv[k][ind_df[k][0]][m1] = log_sum_exp(sIgv, M * M);
			}

			for (int m2 = 0; m2 < M; m2++)
			{
				double sIgv[M * M] = { 0 };
				for (int m1 = 0; m1 < M; m1++)
				{
					for (int m3 = 0; m3 < M; m3++)
					{
						sIgv[m1 * M + m3] = f[1][k][m1][m2][m3] + Ivg[k][ind_df[k][0]][m1] + Ivg[k][ind_df[k][2]][m3];
					}
				}
				Igv[k][ind_df[k][1]][m2] = log_sum_exp(sIgv, M * M);
			}

			for (int m3 = 0; m3 < M; m3++)
			{
				double sIgv[M * M] = { 0 };
				for (int m1 = 0; m1 < M; m1++)
				{
					for (int m2 = 0; m2 < M; m2++)
					{
						sIgv[m1 * M + m2] = f[1][k][m1][m2][m3] + Ivg[k][ind_df[k][0]][m1] + Ivg[k][ind_df[k][1]][m2];
					}
				}
				Igv[k][ind_df[k][2]][m3] = log_sum_exp(sIgv, M * M);
			}
		}

		// Ivg update

		for (int v = 0; v < V; v++)
		{
			double sum0 = 0;
			double sum1 = 0;

			for (int i = 0; i < M; i++)
			{
				sum0 += exp(Igv[ind_dv[v][0]][v][i]);
				sum1 += exp(Igv[ind_dv[v][1]][v][i]);
			}

			sum0 = log(sum0);
			sum1 = log(sum1);

			for (int m = 0; m < M; m++)
			{
				Ivg[ind_dv[v][0]][v][m] = Igv[ind_dv[v][1]][v][m] - sum1;
				Ivg[ind_dv[v][1]][v][m] = Igv[ind_dv[v][0]][v][m] - sum0;
			}
		}
	}

	// Step 3: LLR calculation//TODO

	Q[M][V] = { 0 };

	for (int v = 0; v < V; v++)
	{
		for (int m = 0; m < M; m++)
		{
			Q[m][v] = Ap + Igv[ind_dv[v][0]][v][m] + Igv[ind_dv[v][1]][v][m];
		}
	}

    for(int nresource = 0; nresource < K; nresource++)
	{	
		appLlr[nresource][0][1] = Q[0][ind_df[nresource][0]] - Q[1][ind_df[nresource][0]];
		appLlr[nresource][1][1] = Q[0][ind_df[nresource][1]] - Q[1][ind_df[nresource][1]];
		appLlr[nresource][2][1] = Q[0][ind_df[nresource][2]] - Q[1][ind_df[nresource][2]];
	}
}

