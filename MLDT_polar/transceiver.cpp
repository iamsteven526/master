#define _USE_MATH_DEFINES
#include <iostream>
#include <random>
#include <vector>
#include "polar.h"
#include "parameters.h"
using namespace std;

void Transmitter(PolarCode &polar, vector<vector<uint8_t>> &data, vector<vector<uint8_t>> &codeword, double *preTx, double **tx, double *txFilter)
{
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		for (int i = 0; i < DATA_LEN; i++)
		{
			data.at(nuser).at(i) = rand() % 2;
			//cout << int(data.at(nuser).at(i)) << " ";
		}
		//cout << endl;
		if(DECODING_TYPE)
			codeword.at(nuser) = polar.encode(data.at(nuser), nuser);
		else
			codeword.at(nuser) = polar.encode_bp(data.at(nuser), nuser);
		
		//system("pause");

		Modulator(codeword.at(nuser), tx[nuser]);
		if (!SYNCHRONOUS)
		{
			random_device seed;
			mt19937 generator(seed());
			uniform_real_distribution<double> drift(-DRIFT_RANGE/2., DRIFT_RANGE/2.);
			UpSampling(tx[nuser], preTx);
			SRRCGeneration(txFilter, drift(generator));
			PulseShaping(preTx, tx[nuser], txFilter);
		}
	}
}

void Modulator(vector<uint8_t> &codeword, double *tx)
{
	for (int i = 0; i < CODE_LEN; i++)
	{
		tx[i] = 1 - 2 * codeword[i];
	}
}

void Detector(PolarCode &polar, vector<vector<uint8_t>> &data, vector<vector<double>> &appLlr, long double &errBlock, long double &errBit, double ** app, double ***chcoef, double *llr_spa)
{
	/*int section = 4;
	
	for (int i = 0; i < section; i++)
	{
		int flipping = 2 * (rand() % 2) - 1;
		//cout << flipping << " ";
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			for (int j = 0; j < CODE_LEN / section; j++)
			{
				appLlr[nuser][CODE_LEN / section *i + j] *= flipping;
			}
		}
	}*/
	//cout << endl;
	vector<vector<uint8_t>> decodedResult(NUM_USER, vector<uint8_t>(DATA_LEN));

	if (JOINT == 0)
	{
		
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			if (DECODING_TYPE == 1)
				decodedResult[nuser] = polar.decode_scl_llr(LIST_SIZE, appLlr[nuser], nuser);
			else
				decodedResult[nuser] = polar.decode_BP_sep(appLlr[nuser], nuser);
		}
	}
	else
	{
		
		if (DECODING_TYPE == 1)
			polar.decode_jpscl_llr(LIST_SIZE, app, decodedResult);
		else
			polar.decode_BP_Joint(app,decodedResult);
	}

	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		bool errFlag = false;

		for (int i = NBC; i < DATA_LEN; i++)
		{
			if (decodedResult[nuser][i] != data[nuser][i])
			{
				errFlag = true;
				errBit++;
			}
		}
		if (errFlag) errBlock++;
	}
	//system("pause");

}