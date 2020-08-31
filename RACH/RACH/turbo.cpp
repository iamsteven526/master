#include <iostream>
#include "parameters.h"
using namespace std;

void TrellisConstruction(int ***trellis)
{
	for (int i = 0; i < 2; i++)
	{
		trellis[i][0][0] = trellis[i][0][1] = i;
		trellis[i][1][0] = trellis[i][1][1] = i ^ 1;
	}
}

void BCJR(int ***trellis, double *refLlr, double *rxLlr, double *preLlr, double *appLlr, double **alpha, double **beta, double ***gamma, int sch)
{
	double tempMetric[2];
	//---------- branch metric ----------
	for (int i = 0; i < 2; i++)
	{
		gamma[0][i][0] = refLlr[sch] * (1 - 2 * trellis[i][0][0]);
		gamma[0][i][1] = refLlr[sch] * (1 - 2 * trellis[i][1][0]);
	}
	for (int i = 1; i < FFT_SEGMENT + 1; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			gamma[i][j][0] = 0.5*preLlr[(i - 1)*FFT_POINT + sch] + rxLlr[(i - 1)*FFT_POINT + sch] * (1 - 2 * trellis[j][0][0]);
			gamma[i][j][1] = -0.5*preLlr[(i - 1)*FFT_POINT + sch] + rxLlr[(i - 1)*FFT_POINT + sch] * (1 - 2 * trellis[j][1][0]);
		}
	}
	//---------- initialization ----------
	alpha[1][0] = gamma[0][0][0];
	alpha[1][1] = gamma[0][0][1];
	beta[FFT_SEGMENT + 1][0] = gamma[0][0][0];
	beta[FFT_SEGMENT + 1][1] = gamma[0][1][1];
	//---------- forward metric ----------
	memset(tempMetric, 1000, sizeof(double)*2);
	for (int i = 1; i < FFT_SEGMENT; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				int nextStage = trellis[j][k][1];
				tempMetric[nextStage] = (tempMetric[nextStage] == 1000) ? (gamma[i][j][k] + alpha[i][j]) : LogMaxFunction(gamma[i][j][k] + alpha[i][j], tempMetric[nextStage]);
			}
		}
		memcpy(alpha[i + 1], tempMetric, sizeof(double)*2);
		memset(tempMetric, 1000, sizeof(double)*2);
	}
	//---------- backward metric ----------
	for (int i = FFT_SEGMENT; i > 1; i--)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				for (int input = 0; input < 2; input++)
				{
					if (trellis[k][input][1] == j)
					{
						tempMetric[k] = (tempMetric[k] == 1000) ? gamma[i][k][input] + beta[i + 1][j] : LogMaxFunction(gamma[i][k][input] + beta[i + 1][j], tempMetric[k]);
					}
				}
			}
		}
		memcpy(beta[i], tempMetric, sizeof(double)*2);
		memset(tempMetric, 1000, sizeof(double)*2);
	}
	//---------- a-posteriori LLRs ----------
	for (int i = 1; i < FFT_SEGMENT + 1; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int input = 0; input < 2; input++)
			{
				int nextStage = trellis[j][input][1];
				tempMetric[input] = (tempMetric[input] == 1000) ? beta[i + 1][nextStage] + gamma[i][j][input] + alpha[i][j] : LogMaxFunction(tempMetric[input], beta[i + 1][nextStage] + gamma[i][j][input] + alpha[i][j]);
			}
		}
		appLlr[(i - 1)*FFT_POINT + sch] = tempMetric[0] - tempMetric[1];
		memset(tempMetric, 1000, sizeof(double)*2);
	}
}

double LogMaxFunction(double x, double y)
{
	return (x > y) ? x + log(1. + exp(-fabs(x - y))) : y + log(1. + exp(-fabs(x - y)));
}

void TurboProcessor(LDPC &ldpc, int ***trellis, double *tempLlr, double **refLlr, double *rxLlr, double *preLlr, double **appLlr, double **alpha, double **beta, double ***gamma)
{
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		memset(preLlr, 0, sizeof(double)*CODE_LEN);
		memcpy(rxLlr, appLlr[nuser], sizeof(double)*CODE_LEN);
		for (int it = 1; it <= TURBO_IT; it++)
		{
			for (int sch = 0; sch < FFT_POINT; sch++)
			{
				BCJR(trellis, refLlr[nuser], rxLlr, preLlr, tempLlr, alpha, beta, gamma, sch);
			}
			ldpc.SPA(tempLlr, appLlr[nuser], LDPC_IT);
			bool flag = false;
			for (int i = 0; i < ldpc.h_row; i++)
			{
				int temp = 0;
				for (int j = 0; j < ldpc.row_w[i]; j++)
				{
					temp ^= HARD(appLlr[ldpc.row_edge[i][j]]);
				}
				if (temp == 1)
				{
					flag = true;
					break;
				}
			}
			if (!flag || it == TURBO_IT) break;
			for (int i = 0; i < CODE_LEN; i++)
			{
				preLlr[i] = appLlr[nuser][i] - tempLlr[i];
			}
		}
	}
}