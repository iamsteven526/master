#include <iostream>
#include "parameters.h"
#include <cmath>
#include <cstring>
using namespace std;

void TrellisConstruction(int ***trellis)
{
	for (int s = 0; s < 2; s++)
	{
		trellis[s][0][0] = trellis[s][0][1] = s;
		trellis[s][1][0] = trellis[s][1][1] = s ^ 1;
	}
}

void BCJR(int ***trellis, double *rxLlr, double *appLlr, double **alpha, double **beta, double ***gamma)
{
	double tempMetric[2];
	//---------- branch metric ----------
	for (int i = 0; i < BLOCK_LEN; i++)
	{
		for (int s = 0; s < 2; s++)
		{
			gamma[i][s][0] = rxLlr[i] * (1 - 2 * trellis[s][0][0]);
			gamma[i][s][1] = rxLlr[i] * (1 - 2 * trellis[s][1][0]);
		}
	}
	//---------- initialization ----------
	alpha[1][0] = gamma[0][0][0];
	alpha[1][1] = gamma[0][0][1];
	beta[BLOCK_LEN][0] = gamma[0][0][0];
	beta[BLOCK_LEN][1] = gamma[0][1][1];
	//---------- forward metric ----------
	memset(tempMetric, 1000, sizeof(double)*2);
	for (int i = 1; i < BLOCK_LEN - 1; i++)
	{
		for (int s = 0; s < 2; s++)
		{
			for (int input = 0; input < 2; input++)
			{
				int nextStage = trellis[s][input][1];
				tempMetric[nextStage] = (tempMetric[nextStage] == 1000) ? (gamma[i][s][input] + alpha[i][s]) : LogMaxFunction(gamma[i][s][input] + alpha[i][s], tempMetric[nextStage]);
			}
		}
		memcpy(alpha[i + 1], tempMetric, sizeof(double)*2);
		memset(tempMetric, 1000, sizeof(double)*2);
	}
	//---------- backward metric ----------
	for (int i = BLOCK_LEN - 1; i > 1; i--)
	{
		for (int s = 0; s < 2; s++)
		{
			for (int j = 0; j < 2; j++)
			{
				for (int input = 0; input < 2; input++)
				{
					if (trellis[j][input][1] == s)
					{
						tempMetric[j] = (tempMetric[j] == 1000) ? (gamma[i][j][input] + beta[i + 1][s]) : LogMaxFunction(gamma[i][j][input] + beta[i + 1][s], tempMetric[j]);
					}
				}
			}
		}
		memcpy(beta[i], tempMetric, sizeof(double)*2);
		memset(tempMetric, 1000, sizeof(double)*2);
	}
	//---------- a-posteriori LLRs ----------
	for (int i = 1; i < BLOCK_LEN; i++)
	{
		for (int s = 0; s < 2; s++)
		{
			for (int input = 0; input < 2; input++)
			{
				int nextStage = trellis[s][input][1];
				tempMetric[input] = (tempMetric[input] == 1000) ? beta[i + 1][nextStage] + gamma[i][s][input] + alpha[i][s] : LogMaxFunction(tempMetric[input], beta[i + 1][nextStage] + gamma[i][s][input] + alpha[i][s]);
			}
		}
		appLlr[i] = tempMetric[0] - tempMetric[1];
		memset(tempMetric, 1000, sizeof(double)*2);
	}
}

double LogMaxFunction(double x, double y)
{
	return (x > y) ? x + log(1. + exp(-fabs(x - y))) : y + log(1. + exp(-fabs(x - y)));
}
