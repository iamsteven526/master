#define _USE_MATH_DEFINES
#include <iostream>
#include "parameters.h"
#include <cstring>
#include <math.h>
#include <stddef.h>
using namespace std;

static inline double max_element(double const *ar, size_t size)
{
    double max = ar[0];

    for (int i = 1; i < size; i++)
    {
        if (ar[i] > max)
        {
            max = ar[i];
        }
    }

    return max;
}

static inline double log_sum_exp(double const *x, size_t size)
{
    double xm = max_element(x, size);

    double sum = 0;
    double x_arg = 0;

    #pragma omp simd reduction(+:sum) private(x_arg)
    for (int i = 0; i < size; i++)
    {
        x_arg = x[i] - xm;
        sum += exp(x_arg);
    }

    double log_sum = xm + log(sum);

    return log_sum;
}
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
			/*
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
			*/
		}
		//----------message passing--------------

		int V = 4;       // number of users (layers)
		int M = 2;       // number of codewords in each codebook
		int K = 4;       // number of orthogonal resources
		//#define B   12      // number of bits in SCMA signal
		int df = 3;       // number of branches arriving to a resource node
		int dv = 3;       // number of branches from user node

		int F[4][4]  = { {1,0,1,1},{0,1,1,1},{1,1,1,0},{1,1,0,1} };
		int ind_df[4][3]   = {0};
		int ind_dv[4][3]   = {0};	
		for (int k = 0; k < K; k++)
		{
			int idx = 0;
			for (int v = 0; v < V; v++)
			{
				if (F[k][v] == 1)
				{
					ind_df[k][idx] = v;
					idx++;
				}
			}
		}

		for (int v = 0; v < V; v++)
		{
			int idx = 0;
			for (int k = 0; k < K; k++)
			{
				if (F[k][v] == 1)
				{
					ind_dv[v][idx] = k;
					idx++;
				}
			}
		}
		#pragma omp parallel for
		for (int n = 0; n < 1; n++)
		{
			// Step 1: Initial calculations

			double f[4][2][2][2] = {0};

			for (int k = 0; k < K; k++)
			{
				for(int m1 = 0; m1 < M; m1++)
				{
					for (int m2 = 0; m2 < M; m2++)
					{
						for (int m3 = 0; m3 < M; m3++)
						{
							if( k == 0 || k == 1){
								f[k][m1][m2][m3]   = log(app[(k/2)][2*m+(k%2)][(4*m1+2*m2+m3)]);
							}
							else{
								f[k][m1][m2][m3]   = log(app[(k/2)][2*m+(k%2)][(4*m3+2*m1+m2)]);
							}
						}
					}
				}
			}


			double Ap           = log(1.0/M);
			double Igv[4][4][2] = {0};
			double Ivg[4][4][2] = {0};

			for (int k = 0; k < K; k++)
			{
				for (int v = 0; v < V; v++)
				{
					for (int mm = 0; mm < M; mm++)
					{
						Ivg[k][v][mm] = Ap;
					}
				}
			}

			// Step 2: Iterative procedure
			int Niter = 10;
			for (int iter = 0; iter < Niter; iter++)
			{
				// Igv update

				for (int k = 0; k < K; k++)
				{

					for (int m1 = 0; m1 < M; m1++)
					{
						double sIgv[2*2] __attribute__((aligned(64))) = {0};
						for (int m2 = 0; m2 < M; m2++)
						{
							for (int m3 = 0; m3 < M; m3++)
							{
								sIgv[m2*M+m3] = f[k][m1][m2][m3] + Ivg[k][ind_df[k][1]][m2] + Ivg[k][ind_df[k][2]][m3];
							}
						}
						Igv[k][ind_df[k][0]][m1] = log_sum_exp(sIgv, M*M);
					}

					for (int m2 = 0; m2 < M; m2++)
					{
						double sIgv[2*2] __attribute__((aligned(64))) = {0};
						for (int m1 = 0; m1 < M; m1++)
						{
							for(int m3 = 0; m3 < M; m3++)
							{
								sIgv[m1*M+m3] = f[k][m1][m2][m3] + Ivg[k][ind_df[k][0]][m1] + Ivg[k][ind_df[k][2]][m3];
							}
						}
						Igv[k][ind_df[k][1]][m2] = log_sum_exp(sIgv, M*M);
					}

					for (int m3 = 0; m3 < M; m3++)
					{
						double sIgv[2*2] __attribute__((aligned(64))) = {0};
						for (int m1 = 0; m1 < M; m1++)
						{
							for (int m2 = 0; m2 < M; m2++)
							{
								sIgv[m1*M+m2] = f[k][m1][m2][m3] + Ivg[k][ind_df[k][0]][m1] + Ivg[k][ind_df[k][1]][m2];
							}
						}
						Igv[k][ind_df[k][2]][m3] = log_sum_exp(sIgv, M*M);
					}
				}

				// Ivg update

				for (int v = 0; v < V; v++)
				{
					double sum0 = 0;
					double sum1 = 0;
					double sum2 = 0;

					for (int i = 0; i < M; i++)
					{
						sum0 += exp(Igv[ind_dv[v][0]][v][i]);
						sum1 += exp(Igv[ind_dv[v][1]][v][i]);
						sum2 += exp(Igv[ind_dv[v][2]][v][i]); //20201106
					}

					sum0 = log(sum0);
					sum1 = log(sum1);
					sum2 = log(sum2); // 20201106

					for (int mm = 0; mm < M; mm++)
					{
						Ivg[ind_dv[v][0]][v][mm] = Igv[ind_dv[v][1]][v][mm] - sum1 + Igv[ind_dv[v][2]][v][mm] - sum2;
						Ivg[ind_dv[v][1]][v][mm] = Igv[ind_dv[v][0]][v][mm] - sum0 + Igv[ind_dv[v][2]][v][mm] - sum2;
						Ivg[ind_dv[v][2]][v][mm] = Igv[ind_dv[v][0]][v][mm] - sum0 + Igv[ind_dv[v][1]][v][mm] - sum1;
					}
				}
			}
			// Step 3: LLR calculation

			double Q[2][4] = {0};

			for (int v = 0; v < V; v++)
			{
				for (int mm = 0; mm < M; mm++)
				{
					Q[mm][v] = Ap + Igv[ind_dv[v][0]][v][mm] + Igv[ind_dv[v][1]][v][mm] + Igv[ind_dv[v][2]][v][mm];
					//cout << m << "   " << v << "   " << Q[m][v] << endl;
				}
			}
			//appLlr[(j + 1) % NUM_USER][1]
			for (int v = 0; v < V; v++)
			{
				appLlr[(v/2)][2*m+(v%2)] = log(exp(Q[0][v])/exp(Q[1][v]));
				//cout << appLlr[(v/2)][(v%2)] << endl;
				//LLR[2*v][n]     = log((exp(Q[0][v]) + exp(Q[1][v]))/((exp(Q[2][v]) + exp(Q[3][v]))));
				//LLR[2*v + 1][n] = log((exp(Q[0][v]) + exp(Q[2][v]))/((exp(Q[1][v]) + exp(Q[3][v]))));
			}

		}




	}
}