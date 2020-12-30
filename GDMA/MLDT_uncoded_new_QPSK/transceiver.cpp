#define _USE_MATH_DEFINES
#include <iostream>
#include <random>
#include "parameters.h"
#include <algorithm>
#include <cstring>
using namespace std;

void Modulator(int *data, double *tx, int known_drift);

void Transmitter(int** data, double* preTx, double** tx, double* txFilter,double *rxFilter,double **pilot,double *prePilot, int *known_drift)
{
	// known_drift is the bonder of each range
	
	if (!COLLISION)
	{
		known_drift[0] = 0;
		for (int i = 1; i < NUM_USER; i++){
			known_drift[i] = rand() % BLOCK_LEN;
		}
		
		sort(known_drift, known_drift + NUM_USER);
		known_drift[NUM_USER] = BLOCK_LEN;
		for (int i = 1; i < NUM_USER; i++)
			known_drift[i + NUM_USER] = known_drift[i] + BLOCK_LEN;
	}
	
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		for (int i = 0; i < BLOCK_LEN; i++)
		{
			data[nuser][i] = rand() % 2;
			//cout << data[nuser][i];
			if (!SYNCHRONOUS)
			{
				for (int j = 0; j < SPREAD_LEN; j++)
				{
					pilot[nuser][i * SPREAD_LEN + j] = rand() % 2;
				}
			}
		}
		if (Scrambler)
			scrambler(data[nuser]);

		Modulator(data[nuser], tx[nuser], known_drift[nuser]);
		
		if (!SYNCHRONOUS)
		{
			random_device seed;
			mt19937 generator(seed());
			uniform_real_distribution<double> uniform(-0.5 * DRIFT_RANGE, 0.5 * DRIFT_RANGE);
			UpSampling(tx[nuser], preTx, pilot[nuser], prePilot);
			SRRCGeneration(txFilter, uniform(generator));
			PulseShaping(preTx, tx[nuser], txFilter, rxFilter, prePilot, pilot[nuser]);
		}
	}
}

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
			if (DIFF_ENC && (DIFF_RX_SCHEME == 1)) BCJR(trellis, appLlr[nuser], appLlr[nuser], alpha, beta, gamma);
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

		
		for (int nuser = 0; nuser < NUM_USER; nuser++)
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