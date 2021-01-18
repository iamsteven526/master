#define _USE_MATH_DEFINES
#include <iostream>
#include <random>
#include "ldpc.h"
#include "polar.h"
#include "parameters.h"
#include <time.h>
#include <unordered_set>
#include <cstring>
#include <math.h>
using namespace std;


void Encoder(LDPC &ldpc, PolarCode &polar,int **data, int **codeword, int **Interleaver)
{
	for (int i = 0; i < CODE_LEN * CODE_AMOUNT; i++)
		for(int j=0;j<NUM_USER;j++)
		Interleaver[j][i] = i;

	
	if (INTERLEAVER_type)
	{
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			srand(time(NULL));
			int reg = 0, p = 0;
			for (int i = 0; i < CODE_LEN; i++) //�H���q1~k�� ���p�өM��i�ӥ洫
			{
				p = (rand() ) % CODE_LEN;
				reg = Interleaver[nuser][i];
				Interleaver[nuser][i] = Interleaver[nuser][p];
				Interleaver[nuser][p] = reg;
			}

			/*for (int i = 0; i < CODE_LEN; i++)
			{
				cout << Interleaver[nuser][i] << " ";
			}
			cout << endl;*/
		}
		//system("pause");


	}
	else
	{
		srand(time(NULL));
		int reg = 0, p = 0;
		

		for (int i = 0, m = 0; i < FFT_POINT; i++)  //FFT_point = 512
		{
			for (int j = 0; j < FFT_SEGMENT; j++)
			{
				//Interleaver[0][m] = FFT_POINT*(m/FFT_POINT) + ((5*m)%FFT_POINT);
				TODO!!!!!
				Interleaver[0][m] = (FFT_POINT)*((m%CODE_LEN)/(FFT_POINT/CODE_AMOUNT)) + ((160*m + (m/CODE_LEN))%FFT_POINT);
				//cout << (FFT_POINT)*((m%CODE_LEN)/(FFT_POINT/CODE_AMOUNT)) + ((160*m + (m/CODE_LEN))%FFT_POINT) << endl;
				m++;
			}
		}
        //interleaver in a FFT_POINT*4(64)
		//srand(time(NULL));
		reg = 0;
		p = 0;
		
		for (int j = 0; j < FFT_SEGMENT/4; j++)
		{
			for (int i = 0; i < FFT_POINT*4; i++) 
			{
				p = (rand() ) % (FFT_POINT)*4;
				// not corrrect
				/*
				reg = Interleaver[0][(FFT_POINT)*4*j + i];
				Interleaver[0][(FFT_POINT)*4*j + i] = Interleaver[0][(FFT_POINT)*4*j + p];
				Interleaver[0][(FFT_POINT)*4*j + p] = reg;
				*/
			    //correct
			    /*
				reg = Interleaver[0][j + FFT_POINT*i];
				Interleaver[0][j + FFT_POINT*i] = Interleaver[0][j + FFT_POINT*p];
				Interleaver[0][j + FFT_POINT*p] = reg;
				*/
				
			}
		}
		
		//end interleaver

		for (int i = 1; i < NUM_USER; i++)
		{
			for (int j = 0; j < CODE_LEN * CODE_AMOUNT; j++)
			{
				//cout << Interleaver[0][j] << endl;
				Interleaver[i][j] = Interleaver[0][j];
			}
		}
	}
	int* cache_data = new int[DATA_LEN];
	int* cache_codeword = new int[CODE_LEN];
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		for(int num_code = 0; num_code < CODE_AMOUNT; num_code++){
			for (int i = 0; i < DATA_LEN; i++)
			{
				data[nuser][num_code*DATA_LEN + i] = rand() % 2;
				cache_data[i] = data[nuser][num_code*DATA_LEN + i];
			}

			if (CH_CODING_TYPE)
			{
				ldpc.Encoder(data[nuser], codeword[nuser]);
				if (DIFF_ENC) DiffEncoding(codeword[nuser]);
			}
			else
			{
				if (POLAR_DECODING_TYPE)
					polar.encode(cache_data, cache_codeword, nuser);
				else
					polar.encode_bp(data[nuser], codeword[nuser], nuser);
			}
			for (int i = 0; i < CODE_LEN; i++)
			{
				codeword[nuser][num_code*CODE_LEN + i] = cache_codeword[i];
			}
		}
		if (INTERLEAVER)
		{
			vector<int> temp_c(CODE_AMOUNT * CODE_LEN);
			for (int i = 0; i < CODE_AMOUNT * CODE_LEN; i++)
			{
				temp_c[i] = codeword[nuser][i];
			}
			
			for (int i = 0; i < CODE_AMOUNT * CODE_LEN; i++)
			{
				codeword[nuser][i] = temp_c[Interleaver[nuser][i]];
			}
		}
	}
	delete[] cache_codeword;
	delete[] cache_data;
}

void DiffEncoding(int *codeword)
{
	for (int i = 0; i < FFT_POINT; i++)
	{
		for (int j = 1; j < FFT_SEGMENT; j++)
		{
			codeword[j*FFT_POINT + i] ^= codeword[(j - 1)*FFT_POINT + i];
		}
	}
}

void Modulator(int **codeword, double ***chip)
{
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		if (CH_CODING_TYPE) // LDPC
		{
			if (DIFF_ENC)
			{
				for (int j = 0; j < FFT_POINT; j++)
				{
					chip[nuser][0][j] = 1.; // reference symbol
				}
			}
			for (int i = DIFF_ENC, m = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
			{
				for (int j = 0; j < FFT_POINT; j++)
				{
					chip[nuser][i][j] = 1. - 2. * codeword[nuser][m++];
				}
			}
		}
		else // PolarCode
		{
		/*	for (int i = 0, m = 0; i < FFT_POINT; i++)
			{
				for (int j = 0; j < FFT_SEGMENT; j++)
				{
					chip[nuser][j][i] = 1. - 2. * codeword[nuser][m++];
				}
			}*/
			for (int i = 0, m = 0; i < FFT_SEGMENT; i++)
			{
				for (int j = 0; j < FFT_POINT; j++)
				{
					chip[nuser][i][j] = 1. - 2. * codeword[nuser][m++];
				}
			}

		}
	}
}

void MultiCarrierMapper(double ***chip, double ****tx)
{
	if (!OVER_SAMPLE)
	{
		double zeros[FFT_POINT] = { 0 }, real[FFT_POINT], imag[FFT_POINT];
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			for (int i = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
			{
				
				FFT(chip[nuser][i], zeros, FFT_POINT, FFT_LAYER, real, imag, 1, 0); // IFFT
				memset(zeros, 0, sizeof(double) * FFT_POINT);
				for (int j = 0; j < FFT_POINT; j++)
				{
					tx[nuser][i][j + CP_LEN][0] = real[j] * sqrt((double)FFT_POINT);
					tx[nuser][i][j + CP_LEN][1] = imag[j] * sqrt((double)FFT_POINT);
					
				}
				for (int j = 0; j < CP_LEN; j++) // CP insertion
				{
					tx[nuser][i][j][0] = tx[nuser][i][FFT_POINT + j][0];
					tx[nuser][i][j][1] = tx[nuser][i][FFT_POINT + j][1];
				}

				if (!CP_TYPE)
				{
					for (int j = 0; j < CS_LEN; j++) // CP insertion
					{
						tx[nuser][i][FFT_POINT + CP_LEN + j][0] = tx[nuser][i][CP_LEN + j][0];
						tx[nuser][i][FFT_POINT + CP_LEN + j][1] = tx[nuser][i][CP_LEN + j][1];
					}
				}
			}
		}
	}
	else
	{
		double zeros[FFT_POINT * OVER_SAMPLE_RATE] = { 0 }, real[FFT_POINT * OVER_SAMPLE_RATE], imag[FFT_POINT * OVER_SAMPLE_RATE];

		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			for (int i = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
			{
				//---Over Sampling
				
				double over_sample[FFT_POINT * OVER_SAMPLE_RATE] = { 0 };
				for (int j = 0; j < FFT_POINT; j++)
				{ 
					over_sample[j] = chip[nuser][i][j];
				}
				//---
				FFT(over_sample, zeros, FFT_POINT * OVER_SAMPLE_RATE, FFT_LAYER + log2(OVER_SAMPLE_RATE), real, imag, 1, 0); // IFFT
				memset(zeros, 0, sizeof(double) * FFT_POINT);
				for (int j = 0; j < FFT_POINT * OVER_SAMPLE_RATE; j++)
				{
					tx[nuser][i][j + CP_LEN * OVER_SAMPLE_RATE][0] = real[j] * sqrt((double)FFT_POINT * OVER_SAMPLE_RATE);
					tx[nuser][i][j + CP_LEN * OVER_SAMPLE_RATE][1] = imag[j] * sqrt((double)FFT_POINT * OVER_SAMPLE_RATE);
				}
				for (int j = 0; j < CP_LEN * OVER_SAMPLE_RATE; j++) // CP insertion
				{
					tx[nuser][i][j][0] = tx[nuser][i][FFT_POINT * OVER_SAMPLE_RATE + j][0];
					tx[nuser][i][j][1] = tx[nuser][i][FFT_POINT * OVER_SAMPLE_RATE + j][1];
				}

				if (!CP_TYPE)
				{
					for (int j = 0; j < CS_LEN * OVER_SAMPLE_RATE; j++) // CS insertion
					{
						tx[nuser][i][FFT_POINT * OVER_SAMPLE_RATE + CP_LEN * OVER_SAMPLE_RATE + j][0] = tx[nuser][i][CP_LEN * OVER_SAMPLE_RATE + j][0];
						tx[nuser][i][FFT_POINT * OVER_SAMPLE_RATE + CP_LEN * OVER_SAMPLE_RATE + j][1] = tx[nuser][i][CP_LEN * OVER_SAMPLE_RATE + j][1];
					}
				}
			}
		}
	}
	
}

void MultiCarrierDemapper(double ***rx, double ***postRx, double *drift)
{
	if (!OVER_SAMPLE)
	{
		double real[FFT_POINT], imag[FFT_POINT];
		for (int i = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
		{
			for (int j = 0; j < FFT_POINT; j++)	// CP removal
			{
				real[j] = rx[i][j + CP_LEN][0];
				imag[j] = rx[i][j + CP_LEN][1];
			}

			FFT(real, imag, FFT_POINT, FFT_LAYER, postRx[i][0], postRx[i][1], 0, 0); // FFT
			for (int j = 0; j < FFT_POINT; j++)
			{
				postRx[i][0][j] /= sqrt((double)FFT_POINT);
				postRx[i][1][j] /= sqrt((double)FFT_POINT);
			}
		}
	}
	else
	{
		double real[FFT_POINT * OVER_SAMPLE_RATE];
		double imag[FFT_POINT * OVER_SAMPLE_RATE];
		
		for (int i = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
		{
			for (int j = 0; j < FFT_POINT * OVER_SAMPLE_RATE; j++)	// CP removal
			{
				real[j] = rx[i][j + CP_LEN * OVER_SAMPLE_RATE][0];
				imag[j] = rx[i][j + CP_LEN * OVER_SAMPLE_RATE][1];
			}

			FFT(real, imag, FFT_POINT * OVER_SAMPLE_RATE, FFT_LAYER + log2(OVER_SAMPLE_RATE), postRx[i][0], postRx[i][1], 0, 0); // FFT

			for (int j = 0; j < FFT_POINT * OVER_SAMPLE_RATE; j++)
			{
				postRx[i][0][j] /= sqrt((double)FFT_POINT * OVER_SAMPLE_RATE);
 				postRx[i][1][j] /= sqrt((double)FFT_POINT * OVER_SAMPLE_RATE);
			}
		}
	}
}

void Detector(LDPC &ldpc, PolarCode &polar, int **data, double **appLlr, double **refLlr, long double *errCount, double** app, int** Interleaver, int* error_bits_count)
{


	if (CH_CODING_TYPE)
	{
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			if (!TURBO_DEC && !JCD)
			{
				if (JOINT_DEC) ldpc.JointSPA(appLlr[nuser], refLlr[nuser], appLlr[nuser], JOINT_IT, INNER_IT);
				else
				{
					if (DIFF_DEC) DiffDecoding(appLlr[nuser], refLlr[nuser]);
					ldpc.SPA(appLlr[nuser], appLlr[nuser], LDPC_IT);
				}
			}
			bool errFlag = false;
			for (int i = 0; i < DATA_LEN; i++)
			{
				if (HARD(appLlr[nuser][i]) != data[nuser][i])
				{
					errFlag = true;
					errCount[0]++;
				}
			}
			if (errFlag) errCount[1]++;
		}
	}
	else
	{
		vector<vector<int>> Interleaver_invert(NUM_USER, vector<int> (CODE_AMOUNT * CODE_LEN));
		vector<double> temp_appLlr(CODE_AMOUNT * CODE_LEN);
		vector<vector<double>> temp_app(CODE_AMOUNT * CODE_LEN, vector<double>(NUM_LEVEL));

		if (INTERLEAVER)
		{
			for (int nuser = 0; nuser < NUM_USER; nuser++)
			{
				for (int i = 0; i < CODE_AMOUNT * CODE_LEN; i++)
				{
					Interleaver_invert[nuser][Interleaver[nuser][i]] = i;
				}
			}
		}

		vector<vector<int>> decodedResult(NUM_USER, vector<int>(CODE_AMOUNT * DATA_LEN));

		double *  cache_appLlr = new double[CODE_LEN];
		vector<int> cache_decodedResult(DATA_LEN);
		if (!JOINT)
		{
			if (INTERLEAVER)
			{
				for (int i = 0; i < NUM_USER; i++)
				{
					for (int j = 0; j < CODE_AMOUNT * CODE_LEN; j++){
						temp_appLlr[j] = appLlr[i][j];
						//cout << temp_appLlr[j] << endl;
					}

					for (int j = 0; j < CODE_AMOUNT * CODE_LEN; j++)
						appLlr[i][j] = temp_appLlr[Interleaver_invert[i][j]];
				}
			}
            //TODO:
			for(int code_num = 0; code_num < CODE_AMOUNT; code_num++){
				for (int nuser = 0; nuser < NUM_USER; nuser++)
				{
					for (int j = 0; j < CODE_LEN; j++){
                        cache_appLlr[j] = appLlr[nuser][code_num * CODE_LEN + j];
				    }


					if (POLAR_DECODING_TYPE == 1)
						cache_decodedResult = polar.decode_scl_llr(LIST_SIZE, cache_appLlr, nuser);
					else
						cache_decodedResult = polar.decode_BP_sep(cache_appLlr, nuser);

					for (int j = 0; j < DATA_LEN; j++){
                        decodedResult[nuser][code_num*DATA_LEN + j] = cache_decodedResult[j];
						//cout << cache_decodedResult[j] << endl;
				    }
				}

			}
		}
		else
		{

			if (POLAR_DECODING_TYPE == 1)
				polar.decode_jpscl_llr(LIST_SIZE, app, decodedResult, Interleaver);
			else
				polar.decode_BP_Joint(app, decodedResult, Interleaver);
		}
		delete[] cache_appLlr;
		
		unordered_set<int> NBC_index;
		if (NBC)
		{
			if (POLAR_DECODING_TYPE)
			{
				NBC_index.insert(0);
				NBC_index.insert(9);
				NBC_index.insert(15);
				NBC_index.insert(20);
				NBC_index.insert(27);
				NBC_index.insert(68);
				NBC_index.insert(77);
				NBC_index.insert(82);
				NBC_index.insert(107);
				NBC_index.insert(110);
				NBC_index.insert(125);
				NBC_index.insert(315);
				NBC_index.insert(323);
				NBC_index.insert(368);
				NBC_index.insert(374);			
			}
			else
			{
				NBC_index.insert(368);
				NBC_index.insert(399);
				NBC_index.insert(336);
				NBC_index.insert(293);
				NBC_index.insert(322);
				NBC_index.insert(495);
				NBC_index.insert(282);
				NBC_index.insert(306);
				NBC_index.insert(415);
				NBC_index.insert(479);
				NBC_index.insert(352);
				NBC_index.insert(431);
				NBC_index.insert(463);
				NBC_index.insert(384);
				NBC_index.insert(447);
				NBC_index.insert(511);
			}
		}
		

		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			bool errFlag = false;
			for (int i = 0; i < CODE_AMOUNT * DATA_LEN; i++)
			{
				//cout << decodedResult[nuser][i] << "    " << int(data[nuser][i]) << endl;
				if (NBC_index.count(i) == 1 && NBC)
					continue;
				//cout << decodedResult[nuser][i] << " " << data[nuser][i] << " - ";
				if (decodedResult[nuser][i] != int(data[nuser][i]))
				{
					//cout << i << endl;
					errFlag = true;
					//error_bits_count[i]++;
					errCount[0]++;
				}
			}
			//cout << endl;
			//system("pause");
			if (errFlag)
			{
				errCount[1]++;
			}
		}
	}
}

void DiffDecoding(double *appLlr, double *refLlr)
{
	for (int i = 0; i < FFT_POINT; i++)
	{
		double ref = refLlr[i];
		for (int j = 0; j < FFT_SEGMENT; j++)
		{
			if (appLlr[j*FFT_POINT + i] * ref > 0)
			{
				ref = appLlr[j*FFT_POINT + i];
				appLlr[j*FFT_POINT + i] = abs(appLlr[j*FFT_POINT + i]);
			}
			else
			{
				ref = appLlr[j*FFT_POINT + i];
				appLlr[j*FFT_POINT + i] = -abs(appLlr[j*FFT_POINT + i]);
			}
		}
	}
}