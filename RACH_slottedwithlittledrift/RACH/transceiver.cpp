#define _USE_MATH_DEFINES
#include <iostream>
#include <random>
#include "ldpc.h"
#include "parameters.h"
#include <math.h>
#include <cstring>
using namespace std;
#pragma warning(disable:4996)
namespace
{
	random_device seed;
	mt19937 generator(seed());
	uniform_real_distribution<double> uniform(0.0, 1.0);
}

void Packet_generater(int *packet_num, int **packet_time, long double& packet_sum, double G, long double* errCount)
{
	
	int effLen = OVER_SAMPLE ? CP_TYPE ? (CP_LEN + FFT_POINT) * OVER_SAMPLE_RATE  : (CP_LEN + CS_LEN + FFT_POINT) * OVER_SAMPLE_RATE : CP_TYPE ? (CP_LEN + FFT_POINT) : (CP_LEN + CS_LEN + FFT_POINT);
	int packet_dur = ((DIFF_ENC + FFT_SEGMENT) * effLen + (PREABLE_LEN + CP_LEN + CS_LEN) * PREABLE) * Unit;
	int frame_dur = TEST_FRAME_SIZE * packet_dur;
	double chip_rate = CHIP_RATE * pow(10, 6);
	
	for (int nuser = 0; nuser < NUM_USER; nuser++)
		packet_num[nuser] = 0;

	if (SLOTED)
	{
		poisson_distribution<int> poisson(G);
		for (int t = 0; t < TEST_FRAME_SIZE; t++)
		{
			int NUM_SLOTED = poisson(generator);
			packet_sum += NUM_SLOTED;
			if (!ALOHA && NUM_SLOTED > 1)
			{
				errCount[1] += NUM_SLOTED;
				continue;
			}
			if (ALOHA)
			{
				if (NUM_SLOTED > NUM_USER)
				{
					errCount[1] += NUM_SLOTED;
					continue;
				}
				else if (!CE_SCHEME && NUM_SLOTED >= 4)
				{
					double x = uniform(generator);
					if ((NUM_SLOTED == 4 && x > 0.0041) || (NUM_SLOTED == 5 && x > 0.4371))
					{
						errCount[1] += NUM_SLOTED;
						continue;
					}
				}
			}
			for (int i = 0; i < NUM_SLOTED; i++)
			{
				packet_time[i][packet_num[i]] = packet_dur * t;
				packet_num[i]++;
			}
		}
	}
	else
	{
		int nuser = 0;
		int last_packet_time = 0;
		exponential_distribution<double> exponential(G);
		poisson_distribution<int> poisson(G);
		vector<int> time_in;
		int time_in_num = 0;
		time_in.push_back(-packet_dur);
        int time_drift = 0,cachetimesum = 0;
		int cache_time_drift = 0;
		for (;;)
		{
			int NUM_SLOTED = poisson(generator);
			//NUM_SLOTED = 2;
			//cout << NUM_SLOTED << " ";
			cache_time_drift = 0;
			for(int p = 0; p < NUM_SLOTED; ++p){
				while(true){
					time_drift = exponential(generator) * packet_dur;
					time_drift = time_drift%(2+1);
					cachetimesum = cache_time_drift + time_drift;
					if(cachetimesum > 2*Unit ){
					//if((cachetimesum % (effLen*Unit)) <= 2*Unit ){
					    time_drift = 0;
					}
					time_drift = 0; //slotted
					cache_time_drift = cache_time_drift + time_drift;
					break;
				}
							
				last_packet_time += time_drift;
				if (last_packet_time + packet_dur > frame_dur)
				{
					break;
				}
				if (ALOHA)
				{
					for (int nuser = 0; nuser < NUM_USER; nuser++)
					{
						if (packet_num[nuser] == 0 || last_packet_time - packet_time[nuser][packet_num[nuser] - 1] >= packet_dur)
						{
							//cout << nuser << " ";
							packet_time[nuser][packet_num[nuser]] = last_packet_time;
							packet_num[nuser]++;
							packet_sum++;
							break;
						}
					}
				}
				else
				{
					time_in.push_back(last_packet_time);
					time_in_num++;
					packet_sum++;
				}
			}
			if (last_packet_time + packet_dur*Unit > frame_dur)
			{
				//cout << endl;//debug
				break;
			}
			last_packet_time += (packet_dur)*Unit;				
		}

		time_in.push_back(frame_dur);
	
		if (!ALOHA)
		{
			for (int i = 1; i <= time_in_num; i++)
			{
				if (time_in[i] - time_in[i - 1] > packet_dur && time_in[i + 1] - time_in[i] > packet_dur)
				{
					packet_time[0][packet_num[0]++] = time_in[i];
				}
				else
					errCount[1]++;
			}
		}
	
		
		// 2 user collision case:
		/*int drift_2 = rand() % ((20*63 + 63 + 8)*Unit);
	
		
		if (drift_2 % (effLen * Unit) < 4 * Unit)
		{
			drift_2 -= drift_2 % (effLen * Unit);
		}
		
		
		packet_num[0] = 1;
		packet_num[1] = 1;
		//packet_num[1] = 0;
		packet_time[0][0] = PREABLE_LEN * PREABLE + 0;
		packet_time[1][0] = PREABLE_LEN * PREABLE + drift_2;
		packet_sum += 2;
		//packet_sum += 1;*/

		//cout << PREABLE_LEN + 0 << " " << PREABLE_LEN + drift_2 << endl;
		
		/*for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			for (int t = 0; t < packet_num[nuser]; t++)
			{
				if (packet_time[nuser][t] % effLen < CP_LEN)
					packet_time[nuser][t] -= packet_time[nuser][t] % effLen;
			}
		}
		*/
	}
	
	
	
	/*for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		if (packet_num[nuser] == 0)
			continue;
		for (int t = 0; t < packet_num[nuser]; t++)
		{
			cout << packet_time[nuser][t] << " ";
		}
		cout << endl;
	}
	system("pause");*/
}

void Encoder(LDPC &ldpc, int ***data, int ***codeword, int *packet_num)
{
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		for (int t = 0; t < packet_num[nuser]; t++)
		{
			for (int i = 0; i < DATA_LEN; i++)
			{
				//data[t][nuser][i] = rand() % 2;
				data[t][nuser][i] = 1;
			}
			ldpc.Encoder(data[t][nuser], codeword[t][nuser]);
			if (DIFF_ENC) DiffEncoding(codeword[t][nuser]);
		}
	}

/*	for (int i = 0; i < CODE_LEN; i++)
	{
		cout << codeword[0][0][i] << " ";
	}
	cout << endl;*/
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

void Modulator(int ***codeword, double ****chip, int *packet_num)
{
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		for (int t = 0; t < packet_num[nuser]; t++)
		{
			if (DIFF_ENC)
			{
				for (int j = 0; j < FFT_POINT; j++)
				{
					chip[t][nuser][0][j] = 1.; // reference symbol
				}
			}
			for (int i = DIFF_ENC, m = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
			{
				for (int j = 0; j < FFT_POINT; j++)
				{
					chip[t][nuser][i][j] = 1. - 2. * codeword[t][nuser][m++];
				}
			}
		}
	}
}

void MultiCarrierMapper(double ****chip, double *****tx, int *packet_num)
{
	if (!OVER_SAMPLE)
	{
		double zeros[FFT_POINT] = { 0 }, real[FFT_POINT], imag[FFT_POINT];
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			for (int t = 0; t < packet_num[nuser]; t++)
			{
				for (int i = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
				{
					FFT(chip[t][nuser][i], zeros, FFT_POINT, FFT_LAYER, real, imag, 1, 0); // IFFT
					memset(zeros, 0, sizeof(double) * FFT_POINT);
					for (int j = 0; j < FFT_POINT; j++)
					{
						tx[t][nuser][i][j + CP_LEN][0] = real[j] * sqrt((double)FFT_POINT);
						tx[t][nuser][i][j + CP_LEN][1] = imag[j] * sqrt((double)FFT_POINT);

					}
					for (int j = 0; j < CP_LEN; j++) // CP insertion
					{
						tx[t][nuser][i][j][0] = tx[t][nuser][i][FFT_POINT + j][0];
						tx[t][nuser][i][j][1] = tx[t][nuser][i][FFT_POINT + j][1];
					}

					if (!CP_TYPE)
					{
						for (int j = 0; j < CS_LEN; j++) // CS insertion
						{
							tx[t][nuser][i][FFT_POINT + CP_LEN + j][0] = tx[t][nuser][i][CP_LEN + j][0];
							tx[t][nuser][i][FFT_POINT + CP_LEN + j][1] = tx[t][nuser][i][CP_LEN + j][1];
						}
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
			for (int t = 0; t < packet_num[nuser]; t++)
			{
				for (int i = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
				{
					//---Over Sampling

					double over_sample[FFT_POINT * OVER_SAMPLE_RATE] = { 0 };
					for (int j = 0; j < FFT_POINT; j++)
					{
						over_sample[j] = chip[t][nuser][i][j];
					}
					//---
					FFT(over_sample, zeros, FFT_POINT * OVER_SAMPLE_RATE, FFT_LAYER + log2(OVER_SAMPLE_RATE), real, imag, 1, 0); // IFFT
					memset(zeros, 0, sizeof(double) * FFT_POINT);
					for (int j = 0; j < FFT_POINT * OVER_SAMPLE_RATE; j++)
					{
						tx[t][nuser][i][j + CP_LEN * OVER_SAMPLE_RATE][0] = real[j] * sqrt((double)FFT_POINT * OVER_SAMPLE_RATE);
						tx[t][nuser][i][j + CP_LEN * OVER_SAMPLE_RATE][1] = imag[j] * sqrt((double)FFT_POINT * OVER_SAMPLE_RATE);
					}
					for (int j = 0; j < CP_LEN * OVER_SAMPLE_RATE; j++) // CP insertion
					{
						tx[t][nuser][i][j][0] = tx[t][nuser][i][FFT_POINT * OVER_SAMPLE_RATE + j][0];
						tx[t][nuser][i][j][1] = tx[t][nuser][i][FFT_POINT * OVER_SAMPLE_RATE + j][1];
					}

					if (!CP_TYPE)
					{
						for (int j = 0; j < CS_LEN * OVER_SAMPLE_RATE; j++) // CS insertion
						{
							tx[t][nuser][i][FFT_POINT * OVER_SAMPLE_RATE + CP_LEN * OVER_SAMPLE_RATE + j][0] = tx[t][nuser][i][CP_LEN * OVER_SAMPLE_RATE + j][0];
							tx[t][nuser][i][FFT_POINT * OVER_SAMPLE_RATE + CP_LEN * OVER_SAMPLE_RATE + j][1] = tx[t][nuser][i][CP_LEN * OVER_SAMPLE_RATE + j][1];
						}
					}
				}
			}
		}
	}
	
}

void ReceivingClassify(double** rx, double***** ClassifyRx, int* packet_num, int** packet_time, int ***Cluster_num, int ****Cluster_gain)
{	
	int effLen = OVER_SAMPLE ? CP_TYPE ? (CP_LEN + FFT_POINT) * OVER_SAMPLE_RATE : (CP_LEN + CS_LEN + FFT_POINT) * OVER_SAMPLE_RATE : CP_TYPE ? (CP_LEN + FFT_POINT) : (CP_LEN + CS_LEN + FFT_POINT);
	//---- Cluster num and Cluster gain derive
	vector<vector<bool>> num(NUM_USER, vector<bool>( ( (FFT_SEGMENT + DIFF_ENC) * effLen + (CP_LEN + CS_LEN + PREABLE_LEN) * PREABLE ) * TEST_FRAME_SIZE * Unit));

	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		for (int t = 0; t < packet_num[nuser]; t++)
		{
			for (int i = 0; i < effLen * (FFT_SEGMENT + DIFF_ENC) * Unit; i++)
			{
				num[nuser][packet_time[nuser][t] + i] = 1;
			}
		}
	}

	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		for (int t = 0; t < packet_num[nuser]; t++)
		{
			for (int i = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
			{
				int gain = 0;
				Cluster_num[t][nuser][i] = 0;
				for (int u = 0; u < NUM_USER; u++)
				{
					Cluster_num[t][nuser][i] += num[u][packet_time[nuser][t] + (CP_LEN + i * effLen) * Unit];
					if (num[u][packet_time[nuser][t] + (CP_LEN + i * effLen) * Unit] == 1)
					{
						Cluster_gain[t][nuser][i][gain++] = u;
					}
				}
			}
		}
	}
	num.clear();


	//---- detach the preamble and data
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		for (int t = 0; t < packet_num[nuser]; t++)
		{
			/*
			for (int j = 0; j < PREABLE_LEN * PREABLE; j++)
			{
				cout << rx[packet_time[nuser][t] + j + CP_LEN][0] << " ";
			}
			cout << endl;
			for (int j = 0; j < PREABLE_LEN * PREABLE; j++)
			{
				cout << rx[packet_time[nuser][t] + j + CP_LEN][1]<<" ";
			}
			system("pause");
			*/

			for (int i = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
			{
				for (int j = 0; j < effLen; j++)
				{
					ClassifyRx[t][nuser][i][j][0] = rx[packet_time[nuser][t] + (j + i * effLen + (PREABLE_LEN + CP_LEN + CS_LEN) * PREABLE) * Unit][0];
					ClassifyRx[t][nuser][i][j][1] = rx[packet_time[nuser][t] + (j + i * effLen + (PREABLE_LEN + CP_LEN + CS_LEN) * PREABLE) * Unit][1];
				}	
			}
		}
	}
	
}

void MultiCarrierDemapper(double *****ClassifyRx, double *****postRx, int *packet_num)
{
	/*
	if (EQULIZER)
	{
		vector<vector<vector<double>>> h_matrix(FFT_POINT ,vector<vector<double>> (2 * FFT_POINT, vector<double> (2)));
		vector<vector<vector<double>>> h_zf_matrix(FFT_POINT, vector<vector<double>>(2 * FFT_POINT, vector<double>(2)));
		for (int i = 0; i < FFT_POINT; i++)
		{
			// h1
			for (int j = 0; j < TAP_NUM; j++)
			{
				h_matrix[i][(FFT_POINT - j + 1) % FFT_POINT][0] = h[0][0][j][0];	// real
				h_matrix[i][(FFT_POINT - j + 1) % FFT_POINT][1] = h[0][0][j][1];	// imag
			}

			// h2
			for (int j = 0; j < TAP_NUM; j++)
			{
				h_matrix[i][(FFT_POINT - j + 1) % FFT_POINT + FFT_POINT][0] = h[1][0][j][0];	// real
				h_matrix[i][(FFT_POINT - j + 1) % FFT_POINT + FFT_POINT][1] = h[1][0][j][1];	// imag
			}
		}
	}
	*/

	if (!OVER_SAMPLE)
	{
		double real[FFT_POINT], imag[FFT_POINT];

		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			for (int t = 0; t < packet_num[nuser]; t++)
			{
				for (int i = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
				{
					for (int j = 0; j < FFT_POINT; j++)	// CP removal
					{
						real[j] = ClassifyRx[t][nuser][i][j + CP_LEN][0];
						imag[j] = ClassifyRx[t][nuser][i][j + CP_LEN][1];
					}

					FFT(real, imag, FFT_POINT, FFT_LAYER, postRx[t][nuser][i][0], postRx[t][nuser][i][1], 0, 0); // FFT
					for (int j = 0; j < FFT_POINT; j++)
					{
						postRx[t][nuser][i][0][j] /= sqrt((double)FFT_POINT);
						postRx[t][nuser][i][1][j] /= sqrt((double)FFT_POINT);
					}
				}
			}
		}
	}
	else
	{
		double real[FFT_POINT * OVER_SAMPLE_RATE];
		double imag[FFT_POINT * OVER_SAMPLE_RATE];
		
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			for (int t = 0; t < packet_num[nuser]; t++)
			{
				for (int i = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
				{
					for (int j = 0; j < FFT_POINT * OVER_SAMPLE_RATE; j++)	// CP removal
					{
						real[j] = ClassifyRx[t][nuser][i][j + CP_LEN * OVER_SAMPLE_RATE][0];
						imag[j] = ClassifyRx[t][nuser][i][j + CP_LEN * OVER_SAMPLE_RATE][1];
					}

					FFT(real, imag, FFT_POINT * OVER_SAMPLE_RATE, FFT_LAYER + log2(OVER_SAMPLE_RATE), postRx[t][nuser][i][0], postRx[t][nuser][i][1], 0, 0); // FFT

					for (int j = 0; j < FFT_POINT * OVER_SAMPLE_RATE; j++)
					{
						postRx[t][nuser][i][0][j] /= sqrt((double)FFT_POINT * OVER_SAMPLE_RATE);
						postRx[t][nuser][i][1][j] /= sqrt((double)FFT_POINT * OVER_SAMPLE_RATE);
					}
				}
			}
		}
	}
}

void Detector(LDPC &ldpc, int ***data, double ***appLlr, double ***refLlr, long double *errCount, int *packet_num, bool **estimate_packet_time, double** segment_num, double** segment_error, int **packet_time)
{
	
	/*for (int i = 0; i < CODE_LEN; i++)
	{
		cout << appLlr[0][0][i] << " ";
	}
	system("pause");*/

	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		for (int t = 0; t < packet_num[nuser]; t++)
		{
			//segment_num[nuser][packet_time[nuser][t]]++;
			//cout << "*";
			if (Time_Estimate && estimate_packet_time[nuser][t] == 0)
			{
				errCount[1]++;
				continue;
			}

			if (!TURBO_DEC && !JCD)
			{
				if (JOINT_DEC) ldpc.JointSPA(appLlr[t][nuser], refLlr[t][nuser], appLlr[t][nuser], JOINT_IT, INNER_IT);
				else
				{
					if (DIFF_DEC) DiffDecoding(appLlr[t][nuser], refLlr[t][nuser]);
					ldpc.SPA(appLlr[t][nuser], appLlr[t][nuser], LDPC_IT);
				}
			}
			bool errFlag = false;
			for (int i = 0; i < DATA_LEN; i++)
			{
				if (HARD(appLlr[t][nuser][i]) != data[t][nuser][i])
				{
					errFlag = true;
					break;
					//errCount[0]++;
				}
			}
			if (errFlag)
			{
				//segment_error[nuser][packet_time[1][t]]++;
				errCount[1]++;
				//cout << "+";
				//system("pause");
			}
		}
	}
	//system("pause");

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

void Time_Num_Estimation(double** rx, bool** estimate_packet_time, double& Time_deviation, double& Time_error, int* packet_num, int** packet_time, int* Time_Offset_Error, double***** estimate, double ***ZC_sequence)
{
	int effLen = OVER_SAMPLE ? CP_TYPE ? (CP_LEN + FFT_POINT) * OVER_SAMPLE_RATE : (CP_LEN + CS_LEN + FFT_POINT) * OVER_SAMPLE_RATE : CP_TYPE ? (CP_LEN + FFT_POINT) : (CP_LEN + CS_LEN + FFT_POINT);
	int packet_dur = ((DIFF_ENC + FFT_SEGMENT) * effLen + (PREABLE_LEN + CP_LEN + CS_LEN) * PREABLE ) * Unit;
	int frame_dur = TEST_FRAME_SIZE * packet_dur;

	if (Time_Estimate)
	{
		//---- coarse time offset
		vector<double> cross_correlation_coef;
		vector<int> coarse_offset;
		vector<int> fine_offset;
		int estimate_num = 0;
		for (int r = 0; r < Root_Num; r++)
		{
			double ZC_mean[2] = { 0 };
			double ZC_var[2] = { 0 };

			for (int j = 0; j < PREABLE_LEN; j++)
			{
				ZC_mean[0] += ZC_sequence[r][j][0];
				ZC_mean[1] += ZC_sequence[r][j][1];
			}
			ZC_mean[0] /= PREABLE_LEN;
			ZC_mean[1] /= PREABLE_LEN;

			for (int j = 0; j < PREABLE_LEN; j++)
			{
				ZC_var[0] += pow(ZC_sequence[0][j][0] - ZC_mean[0], 2);
				ZC_var[1] += pow(ZC_sequence[0][j][1] - ZC_mean[1], 2);
			}

			for (int i = 0; i < frame_dur - PREABLE_LEN * Unit; i++)
			{
				double temp[4] = { 0 };
				double means[4] = { 0 };
				//---- cross correlation = conv(rx,filplr(conj(ZC_sequence)))
				for (int j = 0; j < PREABLE_LEN; j++)
				{
					means[2] += rx[i + j * Unit][0];
					means[3] += rx[i + j * Unit][1];
				}
				means[2] /= PREABLE_LEN;
				means[3] /= PREABLE_LEN;

				for (int j = 0; j < PREABLE_LEN; j++)
				{
					temp[0] += (rx[i + j * Unit][0] - means[2]) * (ZC_sequence[r][j][0] - ZC_mean[0]) + (rx[i + j * Unit][1] - means[3]) * (ZC_sequence[r][j][1] - ZC_mean[1]);
					temp[1] += -(rx[i + j * Unit][0] - means[2]) * (ZC_sequence[r][j][1] - ZC_mean[1]) + (rx[i + j * Unit][1] - means[3]) * (ZC_sequence[r][j][0] - ZC_mean[0]);
				}

				for (int j = 0; j < PREABLE_LEN; j++)
				{
					temp[2] += pow(rx[i + j * Unit][0] - means[2], 2);
					temp[3] += pow(rx[i + j * Unit][1] - means[3], 2);
				}

				double t = (pow(temp[0], 2) + pow(temp[1], 2)) / ((temp[2] + temp[3]) * (ZC_var[0] + ZC_var[1])) < NUMERIC_LIMIT ? 0 : (pow(temp[0], 2) + pow(temp[1], 2)) / ((temp[2] + temp[3]) * (ZC_var[0] + ZC_var[1]));
				cross_correlation_coef.push_back(t);
			}

			for (int i = 1; i < frame_dur - PREABLE_LEN * Unit; i++)
			{
				if (cross_correlation_coef[i] > cross_correlation_coef[i - 1] && cross_correlation_coef[i] > cross_correlation_coef[i + 1] && cross_correlation_coef[i] > Threshold)
				{
					coarse_offset.push_back(i - CP_LEN * Unit);
					estimate_num++;
				}
			}
			cross_correlation_coef.clear();
		}

		/*for (int i = 0; i < coarse_offset.size(); i++)
		{
			cout << coarse_offset[i] << " ";
		}
		system("pause");*/


		//---- Fine Time Offset

		vector<vector<double>> estimate_preamble_gain(PREABLE_LEN, vector<double>(2));
		vector<vector<double>> rx_preamble(2, vector<double>(PREABLE_LEN));
		vector<vector<double>> freq_rx_preamble(2, vector<double>(PREABLE_LEN));
		vector<vector<double>> pre_fft_preamble(2, vector<double>(PREABLE_LEN));
		vector<vector<double>> fft_preamble(2, vector<double>(PREABLE_LEN));
		vector<vector<double>> estimate_tap(PREABLE_LEN, vector<double>(2));
		vector<double> tap_e;

		for (int i = 0; i < PREABLE_LEN; i++)
		{
			pre_fft_preamble[0][i] = ZC_sequence[0][i][0];
			pre_fft_preamble[1][i] = ZC_sequence[0][i][1];
		}

		for (int i = 0; i < PREABLE_LEN; i++)
		{
			fft_preamble[0][i] = 0;
			fft_preamble[1][i] = 0;
			for (int j = 0; j < PREABLE_LEN; j++)
			{
				fft_preamble[0][i] += pre_fft_preamble[0][j] * cos(-2 * M_PI * i * j / PREABLE_LEN) - pre_fft_preamble[1][j] * sin(-2 * M_PI * i * j / PREABLE_LEN);
				fft_preamble[1][i] += pre_fft_preamble[1][j] * cos(-2 * M_PI * i * j / PREABLE_LEN) + pre_fft_preamble[0][j] * sin(-2 * M_PI * i * j / PREABLE_LEN);
			}
		}
		for (int t = 0; t < coarse_offset.size(); t++)
		{
			//cout << coarse_offset[t] << " -> ";
			for (int j = 0; j < PREABLE_LEN; j++)
			{
				rx_preamble[0][j] = rx[coarse_offset[t] + j + CP_LEN][0];
				rx_preamble[1][j] = rx[coarse_offset[t] + j + CP_LEN][1];
			}

			//---- DFT of Zadoff Chu
			for (int i = 0; i < PREABLE_LEN; i++)
			{
				freq_rx_preamble[0][i] = 0;
				freq_rx_preamble[1][i] = 0;
				for (int j = 0; j < PREABLE_LEN; j++)
				{
					freq_rx_preamble[0][i] += rx_preamble[0][j] * cos(-2 * M_PI * i * j / PREABLE_LEN) - rx_preamble[1][j] * sin(-2 * M_PI * i * j / PREABLE_LEN);
					freq_rx_preamble[1][i] += rx_preamble[1][j] * cos(-2 * M_PI * i * j / PREABLE_LEN) + rx_preamble[0][j] * sin(-2 * M_PI * i * j / PREABLE_LEN);
				}
			}

			for (int j = 0; j < PREABLE_LEN; j++)
			{
				estimate_preamble_gain[j][0] = (freq_rx_preamble[0][j] * fft_preamble[0][j] + freq_rx_preamble[1][j] * fft_preamble[1][j]) / (pow(fft_preamble[0][j], 2) + pow(fft_preamble[1][j], 2));
				estimate_preamble_gain[j][1] = (freq_rx_preamble[1][j] * fft_preamble[0][j] - freq_rx_preamble[0][j] * fft_preamble[1][j]) / (pow(fft_preamble[0][j], 2) + pow(fft_preamble[1][j], 2));
			}




			//---- IDFT
			double max_tap = 0;
			for (int k = 0; k < PREABLE_LEN; k++)
			{
				if (k >= TAP_NUM + CP_LEN && k < PREABLE_LEN - CS_LEN)
					continue;

				for (int j = 0; j < PREABLE_LEN; j++)
				{
					estimate_tap[k][0] += estimate_preamble_gain[j][0] * cos(2 * M_PI * k * j / PREABLE_LEN) - estimate_preamble_gain[j][1] * sin(2 * M_PI * k * j / PREABLE_LEN);
					estimate_tap[k][1] += estimate_preamble_gain[j][1] * cos(2 * M_PI * k * j / PREABLE_LEN) + estimate_preamble_gain[j][0] * sin(2 * M_PI * k * j / PREABLE_LEN);
				}
				estimate_tap[k][0] /= PREABLE_LEN;
				estimate_tap[k][1] /= PREABLE_LEN;

				tap_e.push_back(pow(estimate_tap[k][0], 2) + pow(estimate_tap[k][1], 2));
				max_tap = pow(estimate_tap[k][0], 2) + pow(estimate_tap[k][1], 2) > max_tap ? pow(estimate_tap[k][0], 2) + pow(estimate_tap[k][1], 2) : max_tap;
			}

			double max_Es = 0;
			double max_l = 0;
			for (int i = 0; i < CS_LEN; i++)
			{
				double temp_Es = 0;
				for (int j = 0; j < TAP_NUM - 1; j++)
				{
					if (tap_e[(-i + j + tap_e.size()) % tap_e.size()] < max_tap * 0.01)
						continue;
					temp_Es += tap_e[(-i + j + tap_e.size()) % tap_e.size()];
				}

				if (temp_Es > max_Es)
				{
					max_Es = temp_Es;
					max_l = i;
				}
			}
			fine_offset.push_back(coarse_offset[t] - max_l);
			tap_e.clear();
		}

		//---- Time_ID_Identify
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			for (int t = 0; t < packet_num[nuser]; t++)
			{
				bool tr = 0;
				int n = 0;
				for (n; n < fine_offset.size(); n++)
				{
					if (0 == fine_offset[n] - packet_time[nuser][t])
					{
						//Time_Offset_Error[20 + estimate_packet_offset[n] - packet_time[nuser][t]] ++;
						tr = 1;
						break;
					}
				}

				if (tr == 1)
				{
					estimate_packet_time[nuser][t] = 1;
					//if(!Phase_Estimate)
					//	Time_deviation += (double)abs(estimate_packet_offset[n] - packet_time[nuser][t]);
				}
				else
				{
					Time_error++;
					estimate_packet_time[nuser][t] = 0;
				}
			}
		}
	}

	if (Phase_Estimate)
	{
		double** rx_preamble = new double* [2];
		double** freq_rx_preamble = new double* [2];
		double** pre_fft_preamble = new double* [2];
		double** fft_preamble = new double* [2];
		double** estimation_gain = new double* [2];

		for (int i = 0; i < 2; i++)
		{
			rx_preamble[i] = new double[PREABLE_LEN];
			freq_rx_preamble[i] = new double[PREABLE_LEN];
			pre_fft_preamble[i] = new double[PREABLE_LEN];
			fft_preamble[i] = new double[PREABLE_LEN];
			estimation_gain[i] = new double[PREABLE_LEN];
		}

		for (int i = 0; i < PREABLE_LEN; i++)
		{
			pre_fft_preamble[0][i] = ZC_sequence[0][i][0];
			pre_fft_preamble[1][i] = ZC_sequence[0][i][1];
		}
		
		

		//---- DFT
		for (int i = 0; i < PREABLE_LEN; i++)
		{
			fft_preamble[0][i] = 0;
			fft_preamble[1][i] = 0;
			for (int j = 0; j < PREABLE_LEN; j++)
			{
				fft_preamble[0][i] += pre_fft_preamble[0][j] * cos(-2 * M_PI * i * j / PREABLE_LEN) - pre_fft_preamble[1][j] * sin(-2 * M_PI * i * j / PREABLE_LEN);
				fft_preamble[1][i] += pre_fft_preamble[1][j] * cos(-2 * M_PI * i * j / PREABLE_LEN) + pre_fft_preamble[0][j] * sin(-2 * M_PI * i * j / PREABLE_LEN);
			}
		}
	

		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			for (int t = 0; t < packet_num[nuser]; t++)
			{
				if (Time_Estimate && estimate_packet_time[nuser][t] == 0)
					continue;

				for (int j = 0; j < PREABLE_LEN; j++)
				{
					rx_preamble[0][j] = rx[packet_time[nuser][t] + j + CP_LEN][0];
					rx_preamble[1][j] = rx[packet_time[nuser][t] + j + CP_LEN][1];
				}
				
				//---- DFT of Zadoff Chu
				for (int i = 0; i < PREABLE_LEN; i++)
				{
					freq_rx_preamble[0][i] = 0;
					freq_rx_preamble[1][i] = 0;
					for (int j = 0; j < PREABLE_LEN; j++)
					{
						freq_rx_preamble[0][i] += rx_preamble[0][j] * cos(-2 * M_PI * i * j / PREABLE_LEN) - rx_preamble[1][j] * sin(-2 * M_PI * i * j / PREABLE_LEN);
						freq_rx_preamble[1][i] += rx_preamble[1][j] * cos(-2 * M_PI * i * j / PREABLE_LEN) + rx_preamble[0][j] * sin(-2 * M_PI * i * j / PREABLE_LEN);
					}
				}

				vector<vector<double>> estimate_preamble_gain(PREABLE_LEN, vector<double>(2));
				vector<vector<double>> estimate_tap(FFT_POINT, vector<double> (2));
		
				for (int j = 0; j < PREABLE_LEN; j++)
				{
					estimate_preamble_gain[j][0] = (freq_rx_preamble[0][j] * fft_preamble[0][j] + freq_rx_preamble[1][j] * fft_preamble[1][j]) / (pow(fft_preamble[0][j], 2) + pow(fft_preamble[1][j], 2));
					estimate_preamble_gain[j][1] = (freq_rx_preamble[1][j] * fft_preamble[0][j] - freq_rx_preamble[0][j] * fft_preamble[1][j]) / (pow(fft_preamble[0][j], 2) + pow(fft_preamble[1][j], 2));
				}
				
				
				//---- IDFT
				for (int k = 0; k < TAP_NUM; k++)
				{
					for (int j = 0; j < PREABLE_LEN; j++)
					{
						estimate_tap[k][0] += estimate_preamble_gain[j][0] * cos(2 * M_PI * k * j / PREABLE_LEN) - estimate_preamble_gain[j][1] * sin(2 * M_PI * k * j / PREABLE_LEN);
						estimate_tap[k][1] += estimate_preamble_gain[j][1] * cos(2 * M_PI * k * j / PREABLE_LEN) + estimate_preamble_gain[j][0] * sin(2 * M_PI * k * j / PREABLE_LEN);
					}
					estimate_tap[k][0] /= PREABLE_LEN;
					estimate_tap[k][1] /= PREABLE_LEN;
				}

				//---- DFT
				for (int i = 0; i < FFT_POINT; i++)
				{
					estimate[t][nuser][0][i][0] = 0;
					estimate[t][nuser][0][i][1] = 0;
					for (int j = 0; j < FFT_POINT; j++)
					{
						estimate[t][nuser][0][i][0] += estimate_tap[j][0] * cos(-2 * M_PI * i * j / FFT_POINT) - estimate_tap[j][1] * sin(-2 * M_PI * i * j / FFT_POINT);
						estimate[t][nuser][0][i][1] += estimate_tap[j][1] * cos(-2 * M_PI * i * j / FFT_POINT) + estimate_tap[j][0] * sin(-2 * M_PI * i * j / FFT_POINT);
					}
				}
				/*
				for (int i = 0; i < FFT_POINT; i++)
				{
					cout << t << " " << estimate[t][nuser][0][i][0] << " ";
				}
				cout << endl;
				for (int i = 0; i < FFT_POINT; i++)
				{
					cout << estimate[t][nuser][0][i][1] << " ";
				}
				
				cout << endl;
				cout << "---------------------------------------------------";
				cout << endl;
				*/
			}
		}
		//system("pause");
		for (int i = 0; i < 2; i++)
		{
			delete[] rx_preamble[i];
			delete[] freq_rx_preamble[i];
			delete[] fft_preamble[i];
			delete[] pre_fft_preamble[i];
			delete[] estimation_gain[i];
		}
		delete[] rx_preamble;
		delete[] freq_rx_preamble;
		delete[] fft_preamble;
		delete[] pre_fft_preamble;
		delete[] estimation_gain;
		
	}
}

void ZC_Generater(double*** ZC_sequence)
{
	int N[3] = { 25, 29, 34 };
	//int N[3] = { 1, 2, 4 };

	for (int r = 0; r < Root_Num; r++)
	{
		for (int i = 0; i < PREABLE_LEN; i++)
		{
			ZC_sequence[r][i][0] = cos(-M_PI * N[r] * i * (i + 1) / PREABLE_LEN);
			ZC_sequence[r][i][1] = sin(-M_PI * N[r] * i * (i + 1) / PREABLE_LEN);
		}
	}
}