#include <iostream>
#include "parameters.h"
#include <cstring>
#include <math.h>
using namespace std;

void Modulator(int *data, double *tx, int known_drift);

void AlamoutiEncoder(int **data, double **tx)
{
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		for(int i=0; i < BLOCK_LEN; ++i){
		//---------- messages generation ----------
		    data[nuser][i] = rand() % 2;
		//---------- Alamouti encoding ----------
		    //tx[2*nuser][i] = 1 - 2 * data[nuser][i];
		}
		for(int i=0; i < BLOCK_LEN; ++i){
		//---------- Alamouti encoding ----------
		    if(i%2 == 0){
				tx[2*nuser][i] = 1 - 2 * data[nuser][i];
                tx[2*nuser+1][i+1] = tx[2*nuser][i];
			}
			else{				
                tx[2*nuser+1][i-1] = 1 - 2 * data[nuser][i]; //-tx[2*nuser][i-1];
				tx[2*nuser][i] = -tx[2*nuser+1][i-1];
			}
		}
	}
}

void SignalCombiner(double **chCoef, double **rx, double ***postRx)
{
	for (int i = 0; i < NUM_USER; i++)
	{
		for(int j = 0; j < BLOCK_LEN/2; ++j){
		//double nFactor = sqrt(pow(chCoef[i][0][0], 2) + pow(chCoef[i][0][1], 2) + pow(chCoef[i][1][0], 2) + pow(chCoef[i][1][1], 2));
		    double nFactor = 1;
			//cout << "chCoef  " << chCoef[0][0] << "  " << chCoef[0][1] << "  " << chCoef[1][0] << "  " << chCoef[1][1] << "  " << chCoef[2][0] << "  " << chCoef[2][1] << "  " << chCoef[3][0] << "  " << chCoef[3][1] << endl;
			//cout << "rx  " << rx[2*j][0] << "  " << rx[2*j][1] << "  " << rx[2*j+1][0] << "  " << rx[2*j][1] << endl;
			postRx[i][2*j][0] = (chCoef[2*i][0] * rx[2*j][0] +  chCoef[2*i][1] * rx[2*j][1] + chCoef[2*i+1][0] * rx[2*j+1][0] +  chCoef[2*i+1][1] * rx[2*j+1][1])  / nFactor;
			postRx[i][2*j][1] = (chCoef[2*i][0] * rx[2*j][1] -  chCoef[2*i][1] * rx[2*j][0] + chCoef[2*i+1][1] * rx[2*j+1][0] -  chCoef[2*i+1][0] * rx[2*j+1][1])  / nFactor;
			postRx[i][2*j+1][0] = (chCoef[2*i+1][0] * rx[2*j][0] + chCoef[2*i+1][1] * rx[2*j][1] - chCoef[2*i][0]*rx[2*j+1][0] - chCoef[2*i][1]*rx[2*j+1][1]) / nFactor;
			postRx[i][2*j+1][1] = (chCoef[2*i+1][0] * rx[2*j][1] - chCoef[2*i+1][1] * rx[2*j][0] - chCoef[2*i][1]*rx[2*j+1][0] + chCoef[2*i][0]*rx[2*j+1][1]) / nFactor;
			//cout << "postRx  " << i << " " << postRx[i][2*j][0] << "  " << postRx[i][2*j][1] << "  " << postRx[i][2*j+1][0] << "  " << postRx[i][2*j+1][1] << endl;

			/*
		    postRx[i][0][0] = (chCoef[i][0][0] * rx[0][0] + chCoef[i][0][1] * rx[0][1] + chCoef[i][1][0] * rx[1][0] + chCoef[i][1][1] * rx[1][1]) / nFactor;
		    postRx[i][0][1] = (chCoef[i][0][0] * rx[0][1] - chCoef[i][0][1] * rx[0][0] + chCoef[i][1][1] * rx[1][0] - chCoef[i][1][0] * rx[1][1]) / nFactor;
		    postRx[i][1][0] = (chCoef[i][1][0] * rx[0][0] + chCoef[i][1][1] * rx[0][1] - chCoef[i][0][0] * rx[1][0] - chCoef[i][0][1] * rx[1][1]) / nFactor;
		    postRx[i][1][1] = (chCoef[i][1][0] * rx[0][1] - chCoef[i][1][1] * rx[0][0] - chCoef[i][0][1] * rx[1][0] + chCoef[i][0][0] * rx[1][1]) / nFactor;
			*/
		    //cout << postRx[i][2*j][0] << "   " << postRx[i][2*j][1] << "   " << postRx[i][2*j+1][0] << "   " << postRx[i][2*j+1][1] << "   " << endl;
		}
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

void SuperLevelSpecification(double **chCoef, double ****supLevel)
{
	if (NUM_USER == 1)
	{
		supLevel[0][0][0][0] = supLevel[0][1][0][0] = pow(chCoef[0][0], 2) + pow(chCoef[0][1], 2) + pow(chCoef[1][0], 2) + pow(chCoef[1][1], 2);
		supLevel[0][0][0][1] = supLevel[0][1][0][1] = 0;
	}
	else if (NUM_USER == 2)
	{
		double temp1[2], temp2[2];
		supLevel[0][0][0][0] = supLevel[0][1][0][0] = pow(chCoef[0][0], 2) + pow(chCoef[0][1], 2) + pow(chCoef[1][0], 2) + pow(chCoef[1][1], 2);
		supLevel[0][0][0][1] = supLevel[0][1][0][1] = 0;
		ComplexConjugate(chCoef[0], temp1);
		ComplexMultiplication(temp1, chCoef[2], temp1);
		ComplexConjugate(chCoef[3], temp2);
		ComplexMultiplication(temp2, chCoef[1], temp2);
		supLevel[0][0][1][0] = temp1[0] + temp2[0];
		supLevel[0][0][1][1] = temp1[1] + temp2[1];
		supLevel[1][1][2][0] = supLevel[0][0][1][0];
		supLevel[1][1][2][1] = supLevel[0][0][1][1];
		ComplexConjugate(supLevel[0][0][1], supLevel[0][1][2]);
		supLevel[1][0][1][0] = supLevel[0][1][2][0];
		supLevel[1][0][1][1] = supLevel[0][1][2][1];
		supLevel[1][0][0][0] = supLevel[1][1][0][0] = pow(chCoef[2][0], 2) + pow(chCoef[2][1], 2) + pow(chCoef[3][0], 2) + pow(chCoef[3][1], 2);
		supLevel[1][0][0][1] = supLevel[1][1][0][1] = 0;
		ComplexConjugate(chCoef[0], temp1);
		ComplexMultiplication(temp1, chCoef[3], temp1);
		ComplexConjugate(chCoef[2], temp2);
		ComplexMultiplication(temp2, chCoef[1], temp2);
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

////TODO UPUPUPUPUPUP

void Modulator(int *data, double *tx, int known_drift)

{
	//cout << "mo" << endl;
	if (COLLISION)
	{
		if (SYNCHRONOUS)
		{
			if (DIFF_ENC) data[0] = 0;
			for (int i = 0; i < BLOCK_LEN; i++)
			{
				tx[i] = 1 - 2 * data[i];
			}
			if (DIFF_ENC)
			{
				for (int i = 1; i < BLOCK_LEN; i++)
				{
					tx[i] *= tx[i - 1];
				}
			}
		}
		else
		{
			int spread_code[SPREAD_LEN] = { 0,0,1,1,0,1,1 };
			int temp[BLOCK_LEN] = { 0 };
			if (DIFF_ENC) temp[0] = 0; // reference symbol
			for (int i = 0; i < BLOCK_LEN; i++)
			{
				temp[i] = 1 - 2 * data[i];
			}
			if (DIFF_ENC)
			{
				for (int i = 1; i < BLOCK_LEN; i++)
				{
					temp[i] *= temp[i - 1];
				}
			}
			//--Spread 
			for (int i = 0; i < BLOCK_LEN; i++)
			{
				for (int j = 0; j < SPREAD_LEN; j++)
				{
					tx[i * SPREAD_LEN + j] = temp[i] * (1 - 2 * spread_code[j]);
				}
			}
		}
	}
	else
	{
		vector<double> temp_tx(BLOCK_LEN);
		if (DIFF_ENC) data[0] = 0;
		for (int i = 0; i < BLOCK_LEN; i++)
		{
			temp_tx[i] = 1 - 2 * data[i];
		}
		if (DIFF_ENC)
		{
			for (int i = 1; i < BLOCK_LEN; i++)
			{
				temp_tx[i] *= temp_tx[i - 1];
			}
		}
		
		for (int i = 0; i < BLOCK_LEN; i++)
		{
			tx[i + known_drift] = temp_tx[i];
		}
	}
}

void Detector(int **data, double **appLlr, int ***trellis, double **alpha, double **beta, double ***gamma, long double &error, int *known_drift)
{
	if (COLLISION)
	{
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			if (DIFF_ENC && (DIFF_RX_SCHEME == 0)) DiffDecoding(appLlr[nuser], 0);
			//if (DIFF_ENC && (DIFF_RX_SCHEME == 1)) BCJR(trellis, appLlr[nuser], appLlr[nuser], alpha, beta, gamma);
			if (Scrambler)
			{
				//vector<int> decoded_data(BLOCK_LEN);
				//decoded_data = descrambler(appLlr[nuser]);
				for (int i = DIFF_ENC; i < BLOCK_LEN; i++)
				{
					error += (data[nuser][i] != HARD(appLlr[nuser][i]));
					//error += (data[nuser][i] != decoded_data[i]);
				}
			}
			else
			{
				for (int i = DIFF_ENC; i < BLOCK_LEN; i++)
				{
					//cout << data[nuser][i] << "  dddd  " << appLlr[nuser][i] << endl;
					error += (data[nuser][i] != HARD(appLlr[nuser][i]));
				}
			}
		}
	}
	else
	{
		// �|���}�o�\��
		if (DIFF_ENC && (DIFF_RX_SCHEME == 1))
		{
			cout << "Parameter Set Error !!";
			cout << "In known drift, BCJR, Detector.";
			system("pause");
		}

		
		for (int nuser = 0; nuser < NUM_USER*NUM_TX; nuser++)
		{
			if (DIFF_ENC && (DIFF_RX_SCHEME == 0)) DiffDecoding(appLlr[nuser], known_drift[nuser]);
			for (int i = 0; i < BLOCK_LEN; i++)
			{
				error += (data[nuser][i] != HARD(appLlr[nuser][i + known_drift[nuser]]));
			}
		}
		

		
	}
}

void DiffDecoding(double *appLlr, int known_drift)
{
	double ref = appLlr[0 + known_drift];
	for (int i = 1; i < BLOCK_LEN; i++)
	{
		if (appLlr[i + known_drift] * ref > 0)
		{
			ref = appLlr[i + known_drift];
			appLlr[i + known_drift] = abs(appLlr[i + known_drift]);
		}
		else
		{
			ref = appLlr[i + known_drift];
			appLlr[i + known_drift] = -abs(appLlr[i + known_drift]);
		}
	}
}

void scrambler(int *data)
{
	int m = 3;
	int s[2] = {2,0};
	int state[3] = {0,0,0};
	int reg = 0;

	for (int i = 0; i < BLOCK_LEN; i++)
	{
		int reg = 0;

		for (int j = 0; j < 2; j++)
		{
			reg += state[s[j]];
		}
		data[i] = (data[i] + reg % 2) % 2;
		
		//--- shift
		for (int j = 4; j > 0; j--)
		{
			state[j] = state[j - 1];
		}
		state[0] = data[i];
	}
	
}

vector<int> descrambler(double* appLlr)
{
	int m = 5;
	int s[4] = { 4,3,1,0 };
	vector<int> temp(BLOCK_LEN + m);
	vector<int> data(BLOCK_LEN);

	for (int i = 0; i < BLOCK_LEN; i++)
	{
		data[i] = HARD(appLlr[i]);
	}

	for (int i = m; i < BLOCK_LEN + m; i++)
	{
		temp[i] = data[i - m];
	}


	for (int i = 0; i < BLOCK_LEN; i++)
	{
		int reg = 0;
		for (int j = 0; j < 4; j++)
		{
			reg += temp[s[j] + i];
		}
		data[i] = (data[i] + reg % 2) % 2;
	}
	
	return data;
}