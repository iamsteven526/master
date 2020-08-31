#define _USE_MATH_DEFINES
#include <iostream>
#include "parameters.h"
using namespace std;

void Transmitter(int **data, double **tx)
{
	for (int nuser = 0; nuser < NUM_USER * SCMA_SOURCE / SCMA_USER_SOURCE; nuser++)
	{
		for (int i = 0; i < BLOCK_LEN*Qm; i++)
		{
			data[nuser][i] = rand() % 2;
		}
		Modulator(data[nuser], tx[nuser]);
	}
}

void Modulator(int *data, double *tx)
{
	double E = sqrt(1. / 2.);
	if (DIFF_ENC)
	{
		double mapping[1 << Qm][Qm] = { { E, E },{ -E, E },{ -E, -E },{ E, -E } };
		int quadrant = 0; // initial quadrant
		tx[0] = tx[1] = E; // reference symbol
		for (int i = 1; i < BLOCK_LEN; i++)
		{
			if (data[Qm * i] == 0 && data[Qm * i + 1] == 1) // phase change: 0 degree
			{
				quadrant = (quadrant + 0) % (1 << Qm);
				tx[Qm * i + 0] = mapping[quadrant][0];
				tx[Qm * i + 1] = mapping[quadrant][1];
			}
			else if (data[Qm * i] == 0 && data[Qm * i + 1] == 0) // phase change: 90 degree
			{
				quadrant = (quadrant + 1) % (1 << Qm);
				tx[Qm * i + 0] = mapping[quadrant][0];
				tx[Qm * i + 1] = mapping[quadrant][1];
			}
			else if (data[Qm * i] == 1 && data[Qm * i + 1] == 0) // phase change: 180 degree
			{
				quadrant = (quadrant + 2) % (1 << Qm);
				tx[Qm * i + 0] = mapping[quadrant][0];
				tx[Qm * i + 1] = mapping[quadrant][1];
			}
			else if (data[Qm * i] == 1 && data[Qm * i + 1] == 1) // phase change: 270 degree
			{
				quadrant = (quadrant + 3) % (1 << Qm);
				tx[Qm * i + 0] = mapping[quadrant][0];
				tx[Qm * i + 1] = mapping[quadrant][1];
			}
		}
	}
	else
	{
		for (int i = 0; i < BLOCK_LEN*Qm; i++)
		{
			tx[i] = E * (1 - 2 * data[i]);
		}
	}
}

void Detecter(int **data, double **appLlr, long double &error)
{
	for (int nuser = 0; nuser < NUM_USER * SCMA_SOURCE / SCMA_USER_SOURCE; nuser++)
	{
		if (DIFF_ENC) DiffDecoding(appLlr[nuser]);
		for (int i = DIFF_ENC*Qm; i < BLOCK_LEN*Qm; i++)
		{
			error += (data[nuser][i] != HARD(appLlr[nuser][i]));
		}
	}
}

void DiffDecoding(double *appLlr)
{
	int mapping[1 << Qm][Qm] = { { 0, 0 },{ 1, 0 },{ 1, 1 },{ 0, 1 } };
	int quadrant;
	for (int i = 0; i < (1 << Qm); i++)
	{
		if (HARD(appLlr[0]) == mapping[i][0] && HARD(appLlr[1]) == mapping[i][1])
		{
			quadrant = i; // reference
			break;
		}
	}
	for (int i = 1; i < BLOCK_LEN; i++)
	{
		for (int j = 0; j < (1 << Qm); j++)
		{
			if (HARD(appLlr[Qm * i]) == mapping[j][0] && HARD(appLlr[Qm * i + 1]) == mapping[j][1])
			{
				int deQuadrant = j;
				if (deQuadrant < quadrant) deQuadrant += (1 << Qm);
				if (deQuadrant - quadrant == 0) // phase change: 0 degree
				{
					appLlr[Qm * i + 0] = +abs(appLlr[Qm * i + 0]);
					appLlr[Qm * i + 1] = -abs(appLlr[Qm * i + 1]);
					quadrant = j;
				}
				else if (deQuadrant - quadrant == 1) // phase change: 90 degree
				{
					appLlr[Qm * i + 0] = +abs(appLlr[Qm * i + 0]);
					appLlr[Qm * i + 1] = +abs(appLlr[Qm * i + 1]);
					quadrant = j;
				}
				else if (deQuadrant - quadrant == 2) // phase change: 180 degree
				{
					appLlr[Qm * i + 0] = -abs(appLlr[Qm * i + 0]);
					appLlr[Qm * i + 1] = +abs(appLlr[Qm * i + 1]);
					quadrant = j;
				}
				else if (deQuadrant - quadrant == 3) // phase change: 270 degree
				{
					appLlr[Qm * i + 0] = -abs(appLlr[Qm * i + 0]);
					appLlr[Qm * i + 1] = -abs(appLlr[Qm * i + 1]);
					quadrant = j;
				}
				break;
			}
		}
	}
}