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
    long double xm = max_element(x, size);

    long double sum = 0;
    long double x_arg = 0;

    #pragma omp simd reduction(+:sum) private(x_arg)
    for (int i = 0; i < size; i++)
    {
        x_arg = x[i] - xm;
        sum += exp(x_arg);
    }
    long double log_sum = xm + log(sum);

    return log_sum;
}

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
				/*
				//TODO:::::
				double prob0, prob1;
				//---------- desired signal ---------- (1)
				prob0 = app[j][i][0] + app[j][i][1] + app[j][i][2] + app[j][i][3] + app[j][i][4] + app[j][i][5] + app[j][i][6] + app[j][i][7] + app[j][i][8] + app[j][i][9] + app[j][i][10] + app[j][i][11] + app[j][i][12] + app[j][i][13] + app[j][i][14] + app[j][i][15];
				prob1 = app[j][i][16] + app[j][i][17] + app[j][i][18] + app[j][i][19] + app[j][i][20] + app[j][i][21] + app[j][i][22] + app[j][i][23] + app[j][i][24] + app[j][i][25] + app[j][i][26] + app[j][i][27] + app[j][i][28] + app[j][i][29] + app[j][i][30] + app[j][i][31];
				
				//cout << "user" << j << "tx " << i << "     0: " << prob0 << "         1: " << prob1 << endl;
				
				if (prob0 <= NUMERIC_LIMIT) appLlr[j][i] += -2*LLR_LIMIT;
				else if (prob1 <= NUMERIC_LIMIT) appLlr[j][i] += 2*LLR_LIMIT;
				else appLlr[j][i] += log(prob0 / prob1);
				//---------- interferences ----------  (4)
				
				prob0 = app[j][i][0] + app[j][i][1] + app[j][i][2] + app[j][i][3] + app[j][i][4] + app[j][i][5] + app[j][i][6] + app[j][i][7] + app[j][i][16] + app[j][i][17] + app[j][i][18] + app[j][i][19] + app[j][i][20] + app[j][i][21] + app[j][i][22] + app[j][i][23];
				prob1 = app[j][i][8] + app[j][i][9] + app[j][i][10] + app[j][i][11] + app[j][i][12] + app[j][i][13] + app[j][i][14] + app[j][i][15] + app[j][i][24] + app[j][i][25] + app[j][i][26] + app[j][i][27] + app[j][i][28] + app[j][i][29] + app[j][i][30] + app[j][i][31];
				if (prob0 <= NUMERIC_LIMIT) appLlr[(3-j)/3][0] += -2*LLR_LIMIT;
				else if (prob1 <= NUMERIC_LIMIT) appLlr[(3-j)/3][0] += 2*LLR_LIMIT;
				else appLlr[(3-j)/3][0] += log(prob0 / prob1);
                
				prob0 = app[j][i][0] + app[j][i][1] + app[j][i][2] + app[j][i][3] + app[j][i][8] + app[j][i][9] + app[j][i][10] + app[j][i][11] + app[j][i][16] + app[j][i][17] + app[j][i][18] + app[j][i][19] + app[j][i][24] + app[j][i][25] + app[j][i][26] + app[j][i][27];
				prob1 = app[j][i][4] + app[j][i][5] + app[j][i][6] + app[j][i][7] + app[j][i][12] + app[j][i][13] + app[j][i][14] + app[j][i][15] + app[j][i][20] + app[j][i][21] + app[j][i][22] + app[j][i][23] + app[j][i][28] + app[j][i][29] + app[j][i][30] + app[j][i][31];
				if (prob0 <= NUMERIC_LIMIT) appLlr[(3-j)/3][1] += -2*LLR_LIMIT;
				else if (prob1 <= NUMERIC_LIMIT) appLlr[(3-j)/3][1] += 2*LLR_LIMIT;
				else appLlr[(3-j)/3][1] += log(prob0 / prob1);				

				prob0 = app[j][i][0] + app[j][i][1] + app[j][i][4] + app[j][i][5] + app[j][i][8] + app[j][i][9] + app[j][i][12] + app[j][i][13] + app[j][i][16] + app[j][i][17] + app[j][i][20] + app[j][i][21] + app[j][i][24] + app[j][i][25] + app[j][i][28] + app[j][i][29];
				prob1 = app[j][i][2] + app[j][i][3] + app[j][i][6] + app[j][i][7] + app[j][i][10] + app[j][i][11] + app[j][i][14] + app[j][i][15] + app[j][i][18] + app[j][i][19] + app[j][i][22] + app[j][i][23] + app[j][i][26] + app[j][i][27] + app[j][i][30] + app[j][i][31];
				if (prob0 <= NUMERIC_LIMIT) appLlr[(7-j)/3][0] += -2*LLR_LIMIT;
				else if (prob1 <= NUMERIC_LIMIT) appLlr[(7-j)/3][0] += 2*LLR_LIMIT;
				else appLlr[(7-j)/3][0] += log(prob0 / prob1);

				prob0 = app[j][i][0] + app[j][i][2] + app[j][i][4] + app[j][i][6] + app[j][i][8] + app[j][i][10] + app[j][i][12] + app[j][i][14] + app[j][i][16] + app[j][i][18] + app[j][i][20] + app[j][i][22] + app[j][i][24] + app[j][i][26] + app[j][i][28] + app[j][i][30];
				prob1 = app[j][i][1] + app[j][i][3] + app[j][i][5] + app[j][i][7] + app[j][i][9] + app[j][i][11] + app[j][i][13] + app[j][i][15] + app[j][i][17] + app[j][i][19] + app[j][i][21] + app[j][i][23] + app[j][i][25] + app[j][i][27] + app[j][i][29] + app[j][i][31];
				if (prob0 <= NUMERIC_LIMIT) appLlr[(7-j)/3][1] += -2*LLR_LIMIT;
				else if (prob1 <= NUMERIC_LIMIT) appLlr[(7-j)/3][1] += 2*LLR_LIMIT;
				else appLlr[(7-j)/3][1] += log(prob0 / prob1);
				*/
			}
		}			
	}
    
	//----------message passing--------------

    int V = 6;       // number of users (layers)
    int M = 2;       // number of codewords in each codebook
    int K = 6;       // number of orthogonal resources
    //#define B   12      // number of bits in SCMA signal
    int df = 5;       // number of branches arriving to a resource node
    int dv = 5;       // number of branches from user node

	int F[6][6]  = { {1,0,1,1,1,1},{0,1,1,1,1,1},{1,1,1,0,1,1},{1,1,0,1,1,1},{1,1,1,1,1,0},{1,1,1,1,0,1} };
    int ind_df[6][5]   = {0};
    int ind_dv[6][5]   = {0};	
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

        double f[6][2][2][2][2][2] = {0};
        for (int k = 0; k < K; k++)
        {
            for(int m1 = 0; m1 < M; m1++)
            {
                for (int m2 = 0; m2 < M; m2++)
                {
                    for (int m3 = 0; m3 < M; m3++)
                    {
						for (int m4 = 0; m4 < M; m4++)
                        {
                            for (int m5 = 0; m5 < M; m5++)
                            {
								if( k == 0 || k == 1){
									f[k][m1][m2][m3][m4][m5]   = log(app[(k/2)][(k%2)][(16*m1+8*m2+4*m3+2*m4+m5)]);
								}
								else if(k == 2 || k == 3){
									f[k][m1][m2][m3][m4][m5]   = log(app[(k/2)][(k%2)][(8*m1+4*m2+16*m3+2*m4+m5)]);
								}
								else{
									f[k][m1][m2][m3][m4][m5]   = log(app[(k/2)][(k%2)][(8*m1+4*m2+2*m3+m4+16*m5)]);
								}
								//cout << f[k][m1][m2][m3][m4][m5] << endl;
							}
						}
                    }
                }
            }
        }


        double Ap           = log(1.0/M);
        long double Igv[6][6][2] = {0};
        long double Ivg[6][6][2] = {0};

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

        // Step 2: Iterative procedure
        int Niter = 10;
        for (int iter = 0; iter < Niter; iter++)
        {
            // Igv update

            for (int k = 0; k < K; k++)
            {
                for (int m1 = 0; m1 < M; m1++)
                {
                    double sIgv[2*2*2*2] __attribute__((aligned(64))) = {0};
                    for (int m2 = 0; m2 < M; m2++)
                    {
                        for (int m3 = 0; m3 < M; m3++)
                        {
                            for (int m4 = 0; m4 < M; m4++)
							{
                                for (int m5 = 0; m5 < M; m5++)
								{
                                    sIgv[m2*M*M*M+m3*M*M+m4*M+m5] = f[k][m1][m2][m3][m4][m5] + Ivg[k][ind_df[k][1]][m2] + Ivg[k][ind_df[k][2]][m3] + Ivg[k][ind_df[k][3]][m4] + Ivg[k][ind_df[k][4]][m5];
								    
								}
							}                         		
                        }
                    }
                    Igv[k][ind_df[k][0]][m1] = log_sum_exp(sIgv, M*M*M*M);
                }

                for (int m2 = 0; m2 < M; m2++)
                {
                    double sIgv[2*2*2*2] __attribute__((aligned(64))) = {0};
                    for (int m1 = 0; m1 < M; m1++)
                    {
                        for(int m3 = 0; m3 < M; m3++)
                        {
                            for (int m4 = 0; m4 < M; m4++)
							{
                                for (int m5 = 0; m5 < M; m5++)
								{
                                    sIgv[m1*M*M*M+m3*M*M+m4*M+m5] = f[k][m1][m2][m3][m4][m5] + Ivg[k][ind_df[k][0]][m1] + Ivg[k][ind_df[k][2]][m3] + Ivg[k][ind_df[k][3]][m4] + Ivg[k][ind_df[k][4]][m5];
								}
							} 
                        }
                    }
                    Igv[k][ind_df[k][1]][m2] = log_sum_exp(sIgv, M*M*M*M);
                }

                for (int m3 = 0; m3 < M; m3++)
                {
                    double sIgv[2*2*2*2] __attribute__((aligned(64))) = {0};
                    for (int m1 = 0; m1 < M; m1++)
                    {
                        for (int m2 = 0; m2 < M; m2++)
                        {
                            for (int m4 = 0; m4 < M; m4++)
							{
                                for (int m5 = 0; m5 < M; m5++)
								{
                                    sIgv[m1*M*M*M+m2*M*M+m4*M+m5] = f[k][m1][m2][m3][m4][m5] + Ivg[k][ind_df[k][0]][m1] + Ivg[k][ind_df[k][1]][m2] + Ivg[k][ind_df[k][3]][m4] + Ivg[k][ind_df[k][4]][m5];
								}
							} 
                        }
                    }
                    Igv[k][ind_df[k][2]][m3] = log_sum_exp(sIgv, M*M*M*M);
                }

                for (int m4 = 0; m4 < M; m4++)
                {
                    double sIgv[2*2*2*2] __attribute__((aligned(64))) = {0};
                    for (int m1 = 0; m1 < M; m1++)
                    {
                        for (int m2 = 0; m2 < M; m2++)
                        {
                            for (int m3 = 0; m3 < M; m3++)
							{
                                for (int m5 = 0; m5 < M; m5++)
								{
                                    sIgv[m1*M*M*M+m2*M*M+m3*M+m5] = f[k][m1][m2][m3][m4][m5] + Ivg[k][ind_df[k][0]][m1] + Ivg[k][ind_df[k][1]][m2] + Ivg[k][ind_df[k][2]][m3] + Ivg[k][ind_df[k][4]][m5];
								}
							} 
                        }
                    }
                    Igv[k][ind_df[k][3]][m4] = log_sum_exp(sIgv, M*M*M*M);
                }

                for (int m5 = 0; m5 < M; m5++)
                {
                    double sIgv[2*2*2*2] __attribute__((aligned(64))) = {0};
                    for (int m1 = 0; m1 < M; m1++)
                    {
                        for (int m2 = 0; m2 < M; m2++)
                        {
                            for (int m3 = 0; m3 < M; m3++)
							{
                                for (int m4 = 0; m4 < M; m4++)
								{
                                    sIgv[m1*M*M*M+m2*M*M+m3*M+m4] = f[k][m1][m2][m3][m4][m5] + Ivg[k][ind_df[k][0]][m1] + Ivg[k][ind_df[k][1]][m2] + Ivg[k][ind_df[k][2]][m3] + Ivg[k][ind_df[k][3]][m4];
								}
							} 
                        }
                    }
                    Igv[k][ind_df[k][4]][m5] = log_sum_exp(sIgv, M*M*M*M);
                }				
            }

            // Ivg update

            for (int v = 0; v < V; v++)
            {
                double sum0 = 0;
                double sum1 = 0;
                double sum2 = 0;
				double sum3 = 0;
				double sum4 = 0;

                for (int i = 0; i < M; i++)
                {
                    sum0 += exp(Igv[ind_dv[v][0]][v][i]);
                    sum1 += exp(Igv[ind_dv[v][1]][v][i]);
                    sum2 += exp(Igv[ind_dv[v][2]][v][i]); //20201106
					sum3 += exp(Igv[ind_dv[v][3]][v][i]);
					sum4 += exp(Igv[ind_dv[v][4]][v][i]);
                }

                sum0 = log(sum0);
                sum1 = log(sum1);
                sum2 = log(sum2); // 20201106
				sum3 = log(sum3);
				sum4 = log(sum4);
				//cout << sum0 << sum1 << sum2 << sum3 << sum4 << endl;

                for (int m = 0; m < M; m++)
                {
                    Ivg[ind_dv[v][0]][v][m] = Igv[ind_dv[v][1]][v][m] - sum1 + Igv[ind_dv[v][2]][v][m] - sum2 + Igv[ind_dv[v][3]][v][m] - sum3 + Igv[ind_dv[v][4]][v][m] - sum4;
                    Ivg[ind_dv[v][1]][v][m] = Igv[ind_dv[v][0]][v][m] - sum0 + Igv[ind_dv[v][2]][v][m] - sum2 + Igv[ind_dv[v][3]][v][m] - sum3 + Igv[ind_dv[v][4]][v][m] - sum4;
                    Ivg[ind_dv[v][2]][v][m] = Igv[ind_dv[v][0]][v][m] - sum0 + Igv[ind_dv[v][1]][v][m] - sum1 + Igv[ind_dv[v][3]][v][m] - sum3 + Igv[ind_dv[v][4]][v][m] - sum4;
					Ivg[ind_dv[v][3]][v][m] = Igv[ind_dv[v][0]][v][m] - sum0 + Igv[ind_dv[v][1]][v][m] - sum1 + Igv[ind_dv[v][2]][v][m] - sum2 + Igv[ind_dv[v][4]][v][m] - sum4;
                    Ivg[ind_dv[v][4]][v][m] = Igv[ind_dv[v][0]][v][m] - sum0 + Igv[ind_dv[v][1]][v][m] - sum1 + Igv[ind_dv[v][2]][v][m] - sum2 + Igv[ind_dv[v][3]][v][m] - sum3;         
				}
            }
        }
        // Step 3: LLR calculation

        long double Q[2][6] = {0};

        for (int v = 0; v < V; v++)
        {
            for (int m = 0; m < M; m++)
            {
                Q[m][v] = Ap + Igv[ind_dv[v][0]][v][m] + Igv[ind_dv[v][1]][v][m] + Igv[ind_dv[v][2]][v][m]  + Igv[ind_dv[v][3]][v][m] + Igv[ind_dv[v][4]][v][m];
				if(!((Q[m][v] <= 10000)&&(Q[m][v] >= -10000)) ){
				    //cout << m << "   " << v << "   " << Q[m][v] << "   " << Igv[ind_dv[v][0]][v][m] << endl;
				}
			}
        }
        //appLlr[(j + 1) % NUM_USER][1]
        for (int v = 0; v < V; v++)
        {
			//appLlr[(v/2)][(v%2)] = log(exp(Q[0][v])/exp(Q[1][v]));
			appLlr[(v/2)][(v%2)] = Q[0][v] - Q[1][v];
			if(!((appLlr[(v/2)][(v%2)] <= 10000)&&(appLlr[(v/2)][(v%2)] >= -10000)) ){
			    //cout << Q[0][v] << "   " << Q[1][v] << "   " << appLlr[(v/2)][(v%2)] << endl;
                //cout << sum0 << endl;
			}
			//cout << appLlr[(v/2)][(v%2)] << endl;
            //LLR[2*v][n]     = log((exp(Q[0][v]) + exp(Q[1][v]))/((exp(Q[2][v]) + exp(Q[3][v]))));
            //LLR[2*v + 1][n] = log((exp(Q[0][v]) + exp(Q[2][v]))/((exp(Q[1][v]) + exp(Q[3][v]))));
        }
	}
	
}