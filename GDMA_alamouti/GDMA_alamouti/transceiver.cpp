#include <iostream>
#include "parameters.h"
#include <cstring>
#include <math.h>
using namespace std;

void AlamoutiEncoder(int **data, double ***tx)
{
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		//---------- messages generation ----------
		data[nuser][0] = rand() % 2;
		data[nuser][1] = rand() % 2;
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
		system("pause");
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