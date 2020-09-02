#define _USE_MATH_DEFINES
#include <iostream>
#include "parameters.h"
#include "scma_log_mpa.h"
#include <omp.h>
using namespace std;

//20200211 add NUM_USER == 3 with CE_METHOD = 1
void fourwayMLDT(double variance, double** chCoef, double*** rx, double** app, double** appLlr, double** estimate, int scma_matrix[SCMA_SOURCE][SCMA_SOURCE * NUM_USER / SCMA_USER_SOURCE])
{
	int count = 0;
	double** chCoef2 = new double* [NUM_USER];
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		chCoef2[nuser] = new double[5];
	}


	for (int i = 0; i < BLOCK_LEN; i++)
	{
		for (int source_id = 0; source_id < SCMA_SOURCE; ++source_id) {
			for (int nuser = 0; nuser < NUM_USER * SCMA_SOURCE / SCMA_USER_SOURCE; nuser++)
			{
				if (scma_matrix[source_id][nuser] == 1) {
					chCoef2[count] = chCoef[nuser];
					++count;
				}
			}
			MLDT(variance, chCoef2, rx[source_id], app[source_id], appLlr, estimate, i);

			count = 0;
		}
		MLDTLLR(i, app, appLlr);
	}
	// modify calculate appLlr by MPA
	
}


void MLDT(double variance, double **chCoef, double **rx, double *app, double** appLlr, double **estimate, int i) //chCoef[userid][modulationlevel]
{
	double chCoefreg = 0, test = 0;
	double probzero = 0 ,probone = 0;
	if (CE_METHOD == 0) {
		for (int i = 0; i < NUM_USER; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				chCoef[i][j] = estimate[i][j];
			}
		}
	}
	//for (int i = 0; i < BLOCK_LEN; i++)
	//{
		//---------- APPs ----------//real
		if (NUM_USER >= 1)//20200211 
		{
			for (int appcounts = 0; appcounts < pow(2, (Qm * NUM_USER)); ++appcounts)
			{
				for (int userid = 0; userid < NUM_USER; ++userid)
				{
					for (int qmlevel = 0; qmlevel < Qm; ++qmlevel)
					{
						test = (Qm * NUM_USER) - userid * 2 - qmlevel;//base for real ex:( - + - + - +)(6 5 4 3 2 1)
						chCoefreg += pow(-1, test + 1) * pow(-1, appcounts / int(pow(2, test - 1))) * chCoef[userid][qmlevel];
					}
				}
				//cout << chCoefreg << endl;//debug
				app[appcounts] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoefreg)), 2.) / (2. * variance));
				//app[appcounts] = -pow((rx[i][0] + sqrt(1. / 2.) * (chCoefreg)), 2.) / (2. * variance);
				chCoefreg = 0;
			}
		}
		normalization(app);
		//---------- APPs ----------//imag
		if (NUM_USER >= 1)//20200211
		{
			for (int appcounts = 0; appcounts < pow(2, (Qm * NUM_USER)); ++appcounts)
			{
				for (int userid = 0; userid < NUM_USER; ++userid)
				{
					for (int qmlevel = 0; qmlevel < Qm; ++qmlevel)
					{
						test = (Qm * NUM_USER) - userid * 2 - (1 - qmlevel);//base for imag ex:( - - - - - -)(5 6 3 4 1 2)
						chCoefreg += (-1) * pow(-1, appcounts / int(pow(2, test - 1))) * chCoef[userid][qmlevel];
					}
				}
				app[appcounts] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoefreg)), 2.) / (2. * variance));
				//app[appcounts] += -pow((rx[i][1] + sqrt(1. / 2.) * (chCoefreg)), 2.) / (2. * variance);
				chCoefreg = 0;
			}
		}
		normalization(app);
		//MLDTLLR(i, app, appLlr);//no scma
	//}
}

void MLDTLLR(int i, double** app, double** appLlr)
{
	
	//const int M = 4;
	//const int V = NUM_USER * SCMA_SOURCE / SCMA_USER_SOURCE;
	//const int K = SCMA_SOURCE;
	double Ap = log(1.0 / M);
	double Igv[K][V][M] = { 0 };
	double Ivg[K][V][M] = { 0 };
	int Niter = 10;
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

	//gogogo initial
	double f[K][M][M][M] = { 0 };
	for (int k = 0; k < K; k++)
	{
		for (int m1 = 0; m1 < M; m1++)
		{
			for (int m2 = 0; m2 < M; m2++)
			{
				for (int m3 = 0; m3 < M; m3++)
				{
					f[k][m1][m2][m3] = log(app[k][16 * m1 + 4 * m2 + m3]);//TODO
				}
			}
		}
	}
	int ind_df[SCMA_SOURCE][NUM_USER] = { {1,2,4}, {0,2,5}, {1,3,5}, {0,3,4} };
	int ind_dv[SCMA_SOURCE * NUM_USER / SCMA_USER_SOURCE][SCMA_USER_SOURCE] = { {1,3}, {0,2}, {0,1}, {2,3}, {0,3}, {1,2} };

	//iteration
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
						sIgv[m2 * M + m3] = f[k][m1][m2][m3] + Ivg[k][ind_df[k][1]][m2] + Ivg[k][ind_df[k][2]][m3];
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
						sIgv[m1 * M + m3] = f[k][m1][m2][m3] + Ivg[k][ind_df[k][0]][m1] + Ivg[k][ind_df[k][2]][m3];
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
						sIgv[m1 * M + m2] = f[k][m1][m2][m3] + Ivg[k][ind_df[k][0]][m1] + Ivg[k][ind_df[k][1]][m2];
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

	//gogogo

	// Step 3: LLR calculation//TODO

	double Q[M][V] = { 0 };

	for (int v = 0; v < V; v++)
	{
		for (int m = 0; m < M; m++)
		{
			Q[m][v] = Ap + Igv[ind_dv[v][0]][v][m] + Igv[ind_dv[v][1]][v][m];
		}
	}

	for (int v = 0; v < V; v++)
	{
		appLlr[v][2*i] = log((exp(Q[0][v]) + exp(Q[1][v])) / ((exp(Q[2][v]) + exp(Q[3][v]))));
		appLlr[v][2*i+1] = log((exp(Q[0][v]) + exp(Q[2][v])) / ((exp(Q[1][v]) + exp(Q[3][v]))));
	}
}
/*
void MLDTLLR(int i, double* app, double** appLlr)//no scma
{
	double test = 0;
	double probzero = 0, probone = 0;
	for (int userid = 0; userid < NUM_USER; ++userid) {
		for (int bitnum = 0; bitnum < 2; ++bitnum) {
			for (int p = 0; p < pow(2, Qm * NUM_USER); ++p) {
				test = (Qm * NUM_USER) - userid * 2 - bitnum;//ex(6 5 4 3 2 1)
				if (pow(-1, p / int(pow(2, test - 1))) == 1) probzero += app[p];
				else probone += app[p];
			}
			if (probzero <= NUMERIC_LIMIT) appLlr[userid][2 * i + bitnum] = -LLR_LIMIT;
			else if (probone <= NUMERIC_LIMIT) appLlr[userid][2 * i + bitnum] = LLR_LIMIT;
			else appLlr[userid][2 * i + bitnum] = log(probzero / probone);
			//cout << probzero << endl;
			probzero = 0;
			probone = 0;
		}
	}
}
*/

void normalization(double* app)
{
	double temp = 0;
	for (int j = 0; j < NUM_LEVEL; j++)
	{
		if (app[j] < NUMERIC_LIMIT) app[j] = NUMERIC_LIMIT;
		temp += app[j];
	}
	for (int j = 0; j < NUM_LEVEL; j++)
	{
		app[j] /= temp;
		if (app[j] < NUMERIC_LIMIT) app[j] = NUMERIC_LIMIT;
	}
}