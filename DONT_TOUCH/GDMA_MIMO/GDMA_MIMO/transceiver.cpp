#include <iostream>
#include "parameters.h"
#include <vector>
#include <cstring>
#include <math.h>
using namespace std;

void AlamoutiEncoder(int **data, double ***tx, int **STBCM)
{
	if (NUM_TX == 2)
	{
		STBCM[0][0] =  1; STBCM[0][1] =  2;
		STBCM[1][0] = -2; STBCM[1][1] =  1;
	}
	else if (NUM_TX == 3)
	{
		STBCM[0][0] = 1; STBCM[0][1] = -2; STBCM[0][2] = -3; STBCM[0][3] = -4; 
		STBCM[1][0] = 2; STBCM[1][1] =  1; STBCM[1][2] =  4; STBCM[1][3] =  3;
		STBCM[2][0] = 3; STBCM[2][1] = -4; STBCM[2][2] =  1; STBCM[2][3] =  2;
	}
	else if (NUM_TX == 4)
	{
		STBCM[0][0] =  1; STBCM[0][1] =  2; STBCM[0][2] =  3; STBCM[0][3] =  4;
		STBCM[1][0] = -2; STBCM[1][1] =  1; STBCM[1][2] = -4; STBCM[1][3] =  3;
		STBCM[2][0] = -3; STBCM[2][1] =  4; STBCM[2][2] =  1; STBCM[2][3] = -2;
		STBCM[3][0] = -4; STBCM[3][1] = -3; STBCM[3][2] =  2; STBCM[3][3] =  1;
	}
	else
	{
		cout << "NUM_TX > 4";
		system("pause");
	}
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		//---------- messages generation ----------
		vector<int> codeword(L);
		for (int l = 0; l < L; l++)
		{
			data[nuser][l] = rand() % 2;
			codeword[l] = 1 - 2 * data[nuser][l];
		}
		//---------- Alamouti encoding ---------
		
		for (int i = 0; i < NUM_TX; i++)
		{
			for (int j = 0; j < L; j++)
			{
				tx[nuser][j][i] = (double)STEP(STBCM[i][j])*codeword[abs(STBCM[i][j]) - 1];
			}
		}
	}
}

void SignalCombiner(double ***chCoef, double **rx, double ***postRx, int **SCM)
{
	if (NUM_TX == 2)
	{
		SCM[0][0] = 1; SCM[0][1] = 2;
		SCM[1][0] = 2; SCM[1][1] = -1;
	}
	else if (NUM_TX == 3)
	{
		SCM[0][0] =  1; SCM[0][1] = 2; SCM[0][2] =  3;
		SCM[1][0] = -2; SCM[1][1] = 1; SCM[1][2] =  4; 
		SCM[2][0] = -3; SCM[2][1] = 4; SCM[2][2] =  1;
		SCM[3][0] = -4; SCM[3][1] = 3; SCM[3][2] = -2;
	}
	else if (NUM_TX == 4)
	{
		SCM[0][0] = 1; SCM[0][1] =  2; SCM[0][2] =  3; SCM[0][3] =  4;
		SCM[1][0] = 2; SCM[1][1] = -1; SCM[1][2] = -4; SCM[1][3] =  3;
		SCM[2][0] = 3; SCM[2][1] =  4; SCM[2][2] = -1; SCM[2][3] = -2;
		SCM[3][0] = 4; SCM[3][1] = -3; SCM[3][2] =  2; SCM[3][3] = -1;
	}
	else
	{
		cout << "NUM_TX > 4";
		system("pause");
	}
	

	for (int i = 0; i < NUM_USER; i++)
	{
	/*	double nFactor = 0;
		for (int j = 0; j < NUM_TX; j++)
			nFactor += pow(chCoef[i][j][0], 2) + pow(chCoef[i][j][1], 2);
		nFactor = sqrt(nFactor);
	*/
	//	nFactor = sqrt(pow(chCoef[i][0][0], 2) + pow(chCoef[i][0][1], 2) + pow(chCoef[i][1][0], 2) + pow(chCoef[i][1][1], 2));
		
		for (int j = 0; j < L; j++)
		{
			postRx[i][j][0] = 0;
			postRx[i][j][1] = 0;
			for (int k = 0; k < NUM_TX; k++)
			{
				postRx[i][j][0] += (double)STEP(SCM[j][k]) * rx[abs(SCM[j][k]) - 1][0] * chCoef[i][k][0] - (double)STEP(SCM[j][k]) * rx[abs(SCM[j][k]) - 1][1] * chCoef[i][k][1];
				postRx[i][j][1] += (double)STEP(SCM[j][k]) * rx[abs(SCM[j][k]) - 1][0] * chCoef[i][k][1] + (double)STEP(SCM[j][k]) * rx[abs(SCM[j][k]) - 1][1] * chCoef[i][k][0];
			}
		//	postRx[i][j][0] /= nFactor;
		//	postRx[i][j][1] /= nFactor;
		}
	}
			
	/* 
	Complex value design         
	| x1 -x2* |
	| x2  x1  |
	*/
	/*
	postRx[i][0][0] = (chCoef[i][0][0] * rx[0][0] + chCoef[i][0][1] * rx[0][1] + chCoef[i][1][0] * rx[1][0] + chCoef[i][1][1] * rx[1][1]) / nFactor;
	postRx[i][0][1] = (chCoef[i][0][0] * rx[0][1] - chCoef[i][0][1] * rx[0][0] + chCoef[i][1][1] * rx[1][0] - chCoef[i][1][0] * rx[1][1]) / nFactor;
	postRx[i][1][0] = (chCoef[i][1][0] * rx[0][0] + chCoef[i][1][1] * rx[0][1] - chCoef[i][0][0] * rx[1][0] - chCoef[i][0][1] * rx[1][1]) / nFactor;
	postRx[i][1][1] = (chCoef[i][1][0] * rx[0][1] - chCoef[i][1][1] * rx[0][0] - chCoef[i][0][1] * rx[1][0] + chCoef[i][0][0] * rx[1][1]) / nFactor;
	*/

		
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

void SuperLevelSpecification(double ***chCoef, double ****supLevel, int **SCM)
{
	if (NUM_USER == 1)
	{
		
		//Complex value matrix design
		/*
		supLevel[0][0][0][0] = supLevel[0][1][0][0] = pow(chCoef[0][0][0], 2) + pow(chCoef[0][0][1], 2) + pow(chCoef[0][1][0], 2) + pow(chCoef[0][1][1], 2);
		supLevel[0][0][0][1] = supLevel[0][1][0][1] = 0;
		*/
		// Real value matrix design
		//cout << chCoef[0][0][0] << " " << chCoef[0][0][1] << " " << chCoef[0][1][0] << " " << chCoef[0][1][1] << endl;
		
		//supLevel[0][0][0][0] = supLevel[0][1][0][0] = pow(chCoef[0][0][0], 2) - pow(chCoef[0][0][1], 2) + pow(chCoef[0][1][0], 2) - pow(chCoef[0][1][1], 2);
		//supLevel[0][0][0][1] = supLevel[0][1][0][1] = 2 * chCoef[0][0][0] * chCoef[0][0][1] + 2 * chCoef[0][1][0] * chCoef[0][1][1];
		
		for (int i = 0; i < L; i++)
		{
			supLevel[0][i][0][0] = 0;
			supLevel[0][i][0][1] = 0;
		}

		for (int j = 0; j < NUM_TX; j++)
		{
			supLevel[0][0][0][0] += pow(chCoef[0][j][0], 2) - pow(chCoef[0][j][1], 2);
			supLevel[0][0][0][1] += 2 * (chCoef[0][j][0] * chCoef[0][j][1]);
		}

		for (int i = 1; i < L; i++)
		{
			supLevel[0][i][0][0] = supLevel[0][0][0][0];
			supLevel[0][i][0][1] = supLevel[0][0][0][1];
		}
		
		//cout << supLevel[0][0][0][0] << " " << supLevel[0][0][0][1] << endl;

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
	else
	{
		printf("\nPARAMETER SETTING IS WRONG\n");
		system("pause");
	}
	/*for (int i = 0; i < NUM_USER; i++)
	{
		double nFactor = sqrt(pow(chCoef[i][0][0], 2) + pow(chCoef[i][0][1], 2) + pow(chCoef[i][1][0], 2) + pow(chCoef[i][1][1], 2));
		for (int j = 0; j < NUM_TX; j++)
		{
			for (int k = 0; k < (NUM_USER*NUM_TX - 1); k++)
			{
				supLevel[i][j][k][0] /= nFactor;
				supLevel[i][j][k][1] /= nFactor;
			}
		}
	}*/
}

void Detector(int **data, double **appLlr, long double &error)
{
	//cout << endl;
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		for (int i = 0; i < NUM_TX; i++)
		{
			//cout << appLlr[nuser][i] << " ";
			error += (data[nuser][i] != HARD(appLlr[nuser][i]));
		}
	}
	//system("pause");
	//cout << endl;
}