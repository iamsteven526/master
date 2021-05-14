#define _USE_MATH_DEFINES
#include <iostream>
#include <random>
#include "parameters.h"
using namespace std;

namespace
{
	random_device seed;
	mt19937 generator(seed());
	uniform_real_distribution<double> uniform(0, 1);
	exponential_distribution<double> exponential(1);
	normal_distribution<double> normal(0, 1);
	uniform_real_distribution<double> uniform_drift(-0.5 * DRIFT_RANGE, 0.5 * DRIFT_RANGE);
}

void EnergyProfile(double *****h, double *****H, int *packet_num, int **packet_time, double ***h_preamble, double ***H_preamble)
{
	//---------- channel coefficients in the time domain ----------

	if (CH_TYPE == 1) // quasic-static
	{
		double nFactor = 0;
		for (int i = 0; i < TAP_NUM; i++)
		{
			nFactor += exp(-i);
		}
		for (int i = 0; i < NUM_USER; i++)
		{
			if (packet_num[i] == 0)
				continue;
			for (int j = 0; j < FFT_POINT; j++)
			{
				if (j < TAP_NUM)
				{
					double power = 0.5*exp(-j) / nFactor;
					double phi = uniform(generator) * 2. * M_PI - M_PI;
					h[0][i][0][j][0] = normal(generator);
					h[0][i][0][j][1] = normal(generator);

					h[0][i][0][j][0] = sqrt(power) * (pow(K + 1, -0.5) * h[0][i][0][j][0] + pow(K, 0.5) / pow(K + 1, 0.5) * cos(phi));
					h[0][i][0][j][1] = sqrt(power) * (pow(K + 1, -0.5) * h[0][i][0][j][1] + pow(K, 0.5) / pow(K + 1, 0.5) * sin(phi));
					//cout << h[i][0][j][0] << " " << h[i][0][j][1] << endl;
				}
				else h[0][i][0][j][0] = h[0][i][0][j][1] = 0;
				for (int k = 1; k < FFT_SEGMENT + DIFF_ENC; k++) // block fading
				{
					h[0][i][k][j][0] = h[0][i][0][j][0];
					h[0][i][k][j][1] = h[0][i][0][j][1];
				}
			}
			//system("pause");
		}
	}
	else if (CH_TYPE == 2) // time-varying
	{
		
		double period = 1. / (double)(CHIP_RATE*pow(10, 6)); // sampling duration
		double fd = (UE_SPEED*0.277778)*(CARRIER_FREQ*pow(10, 9)) / (3. * pow(10, 8)); // maximun Doppler frequency shift
		const int M = 16;
		double nFactor = 0;
		for (int i = 0; i < TAP_NUM; i++)
		{
			nFactor += exp(-i);
		}
		for (int i = 0; i < NUM_USER; i++)
		{
			for (int j = 0; j < FFT_POINT; j++)
			{
				if (j < TAP_NUM)
				{
					h_preamble[i][j][0] = 0;
					h_preamble[i][j][1] = 0;

					double power = 0.5*exp(-j) / nFactor;
					double theta = uniform(generator) * 2. * M_PI - M_PI;
					for (int k = 0; k < FFT_SEGMENT + DIFF_ENC; k++)
					{
						h[0][i][k][j][0] = h[0][i][k][j][1] = 0;
					}
					// Derive c_i and c_q
					
					for (int n = 1; n <= M; n++)
					{
						double phi_1 = uniform(generator) * 2. * M_PI - M_PI;
						double phi_2 = uniform(generator) * 2. * M_PI - M_PI;
						double alpha = (2.*M_PI*n - M_PI + theta) / (4.*M);

						h_preamble[i][j][0] += cos(2. * M_PI * fd * 0 * cos(alpha) + phi_1);
						h_preamble[i][j][1] += cos(2. * M_PI * fd * 0 * sin(alpha) + phi_2);

						for (int k = 0; k < FFT_SEGMENT + DIFF_ENC; k++)
						{
							double t = (k * FFT_POINT + PREABLE_LEN) * period;
							h[0][i][k][j][0] += cos(2.*M_PI*fd*t*cos(alpha) + phi_1);
							h[0][i][k][j][1] += cos(2.*M_PI*fd*t*sin(alpha) + phi_2);
						}
					}

					// Derive different channel tap <---racian channel

					double phi = uniform(generator) * 2. * M_PI - M_PI;
					double thi = uniform(generator) * 2. * M_PI - M_PI;
					for (int k = 0; k < FFT_SEGMENT + DIFF_ENC; k++)
					{
						double t = k * FFT_POINT * period;
						h[0][i][k][j][0] = sqrt(power) * (pow(K + 1, -0.5) * h[0][i][k][j][0] * sqrt(2. / M) + pow(K, 0.5) / pow(K + 1, 0.5) * cos(2. * M_PI * fd * t * cos(thi) + phi));
						h[0][i][k][j][1] = sqrt(power) * (pow(K + 1, -0.5) * h[0][i][k][j][1] * sqrt(2. / M) + pow(K, 0.5) / pow(K + 1, 0.5) * sin(2. * M_PI * fd * t * cos(thi) + phi));
					}

				}
				else
				{
					// Derive different channel tap
					for (int k = 0; k < FFT_SEGMENT + DIFF_ENC; k++)
					{
						h[0][i][k][j][0] = h[0][i][k][j][1] = 0;
					}
				}
			}
		}
	}

	//---------- channel coefficients in the frequency domain ----------
	double real[FFT_POINT], imag[FFT_POINT];
	double real_p[FFT_POINT], imag_p[FFT_POINT];
	for (int i = 0; i < NUM_USER; i++)
	{
		for (int j = 0; j < FFT_SEGMENT + DIFF_ENC; j++)
		{
			for (int k = 0; k < FFT_POINT; k++)
			{
				real[k] = h[0][i][j][(k+FFT_POINT-(packet_time[i][0])%(CP_LEN+CS_LEN+FFT_POINT))%FFT_POINT][0];
				imag[k] = h[0][i][j][(k+FFT_POINT-(packet_time[i][0])%(CP_LEN+CS_LEN+FFT_POINT))%FFT_POINT][1];
			}
			FFT(real, imag, FFT_POINT, FFT_LAYER, H[0][i][j][0], H[0][i][j][1], 0, 0);
		}

		/*if (CH_TYPE == 2)
		{
			for (int k = 0; k < PREABLE_LEN; k++)
			{
				real_p[k] = h_preamble[i][k][0];
				imag_p[k] = h_preamble[i][k][1];
			}
			//---- DFT of Zadoff Chu
			for (int k= 0; k < PREABLE_LEN; k++)
			{
				H_preamble[i][0][k] = 0;
				H_preamble[i][1][k] = 0;
				for (int j = 0; j < PREABLE_LEN; j++)
				{
					H_preamble[i][0][k] += real_p[j] * cos(-2 * M_PI * k * j / PREABLE_LEN) - imag_p[j] * sin(-2 * M_PI * k * j / PREABLE_LEN);
					H_preamble[i][1][k] += imag_p[j] * cos(-2 * M_PI * k * j / PREABLE_LEN) + real_p[j] * sin(-2 * M_PI * k * j / PREABLE_LEN);
				}
			}
		}*/
	}

	/*if (CH_TYPE == 2 && !SLIDING) // average the channel coefficients varing with time for comparison 
	{
		for (int i = 0; i < NUM_USER; i++)
		{
			for (int k = 0; k < FFT_POINT; k++)
			{
				for (int j = 1; j < FFT_SEGMENT + DIFF_ENC; j++)
				{
				
					H[0][i][0][0][k] += H[0][i][j][0][k];
					H[0][i][0][1][k] += H[0][i][j][1][k];
				}
			}
			for (int k = 0; k < FFT_POINT; k++)
			{
				H[0][i][0][0][k] /= (double)(FFT_SEGMENT + DIFF_ENC);
				H[0][i][0][1][k] /= (double)(FFT_SEGMENT + DIFF_ENC);
			}
		}
	}*/

	/*
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		if (packet_num[nuser] == 0)
			continue;
		for (int i = 0; i < FFT_POINT; i++)
		{
			cout << H[0][nuser][0][0][i] << " ";
		}
		cout << endl;
		for (int i = 0; i < FFT_POINT; i++)
		{
			cout << H[0][nuser][0][1][i] << " ";
		}
		cout << endl;
		cout << "---------------------------------------------------";
		cout << endl;
	}
	system("pause");
	*/
}


void MultipleAccessChannel(double stdDev, double***** h, double***** tx, double** rx, double* drift, double***** H, int* packet_num, int** packet_time, double* txFilter, double***ZC_sequence, double*** h_preamble)
{
	int effLen = OVER_SAMPLE ? CP_TYPE ? (CP_LEN + FFT_POINT) * OVER_SAMPLE_RATE : (CP_LEN + CS_LEN + FFT_POINT) * OVER_SAMPLE_RATE : CP_TYPE ? (CP_LEN + FFT_POINT) : (CP_LEN + CS_LEN + FFT_POINT);
	int preamble_len = 0;//PREABLE * (PREABLE_LEN + CP_LEN + CS_LEN);
	int packet_dur = ((FFT_SEGMENT + DIFF_ENC) * effLen + preamble_len) * Unit;
	int trunLen = TRUNCATION * UP_RATE;
	vector<vector<vector<double>>> packet(NUM_USER, vector<vector<double>>(packet_dur* TEST_FRAME_SIZE, vector<double>(2)));
	vector<vector<double>> pre_Filter_rx(packet_dur * TEST_FRAME_SIZE, vector<double>(2));
	

	// MAC
	for (int i = 0; i < packet_dur * TEST_FRAME_SIZE; i++)
	{
		rx[i][0] = stdDev * normal(generator); // real
		rx[i][1] = stdDev * normal(generator); // imaginary
		//rx[i][0] = 0;
		//rx[i][1] = 0;
	}


	// Packet Gernatator
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		for (int t = 0; t < packet_num[nuser]; t++)
		{
			vector<vector<double>> pre_packet(packet_dur, vector <double>(2));
			for (int i = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
			{
				for (int j = 0; j < effLen; j++)
				{
					pre_packet[(j + i * effLen) * Unit][0] = tx[t][nuser][i][j][0];
					pre_packet[(j + i * effLen) * Unit][1] = tx[t][nuser][i][j][1];
				}
			}

			if (CH_TYPE == 1)
			{
				for (int i = 0; i < packet_dur; i++)
				{
					if (PULSE_SHAPE)
					{
						for (int j = -trunLen; j < trunLen + UP_RATE; j++)
						{
							packet[nuser][packet_time[nuser][t] + i][0] += txFilter[j + trunLen] * pre_packet[(packet_dur + (i - j)) % packet_dur][0]; // real
							packet[nuser][packet_time[nuser][t] + i][1] += txFilter[j + trunLen] * pre_packet[(packet_dur + (i - j)) % packet_dur][1]; // imaginary
						}
						//	cout << packet[nuser][packet_time[nuser][t] + i][0] << " ";
					}
					else
					{
						packet[nuser][packet_time[nuser][t] + i][0] += pre_packet[i][0]; // real
						packet[nuser][packet_time[nuser][t] + i][1] += pre_packet[i][1]; // imaginary
					}
				}
			}
			else if (CH_TYPE == 2)
			{
				//----fast fading----
				for (int i = 0; i < preamble_len; i++)
				{
					if (i < TAP_NUM * Unit)
					{
						for (int j = 0; j < i / Unit; j++)
						{
							/*pre_Filter_rx[packet_time[nuser][t] + i][0] += pre_packet[i - j * Unit][0] * h_preamble[nuser][j][0];
							pre_Filter_rx[packet_time[nuser][t] + i][0] -= pre_packet[i - j * Unit][1] * h_preamble[nuser][j][1];
							pre_Filter_rx[packet_time[nuser][t] + i][1] += pre_packet[i - j * Unit][0] * h_preamble[nuser][j][1];
							pre_Filter_rx[packet_time[nuser][t] + i][1] += pre_packet[i - j * Unit][1] * h_preamble[nuser][j][0];*/
							pre_Filter_rx[packet_time[nuser][t] + i][0] += pre_packet[i - j * Unit][0] * h[0][nuser][0][j][0];
							pre_Filter_rx[packet_time[nuser][t] + i][0] -= pre_packet[i - j * Unit][1] * h[0][nuser][0][j][1];
							pre_Filter_rx[packet_time[nuser][t] + i][1] += pre_packet[i - j * Unit][0] * h[0][nuser][0][j][1];
							pre_Filter_rx[packet_time[nuser][t] + i][1] += pre_packet[i - j * Unit][1] * h[0][nuser][0][j][0];
						}
					}
					else
					{
						for (int j = 0; j < TAP_NUM; j++)
						{
							/*pre_Filter_rx[packet_time[nuser][t] + i][0] += pre_packet[(i - j) * Unit][0] * h_preamble[nuser][j][0];
							pre_Filter_rx[packet_time[nuser][t] + i][0] -= pre_packet[(i - j) * Unit][1] * h_preamble[nuser][j][1];
							pre_Filter_rx[packet_time[nuser][t] + i][1] += pre_packet[(i - j) * Unit][0] * h_preamble[nuser][j][1];
							pre_Filter_rx[packet_time[nuser][t] + i][1] += pre_packet[(i - j) * Unit][1] * h_preamble[nuser][j][0];*/
							pre_Filter_rx[packet_time[nuser][t] + i][0] += pre_packet[(i - j) * Unit][0] * h[0][nuser][0][j][0];
							pre_Filter_rx[packet_time[nuser][t] + i][0] -= pre_packet[(i - j) * Unit][1] * h[0][nuser][0][j][1];
							pre_Filter_rx[packet_time[nuser][t] + i][1] += pre_packet[(i - j) * Unit][0] * h[0][nuser][0][j][1];
							pre_Filter_rx[packet_time[nuser][t] + i][1] += pre_packet[(i - j) * Unit][1] * h[0][nuser][0][j][0];
						}
					}
				}

				for (int i = 0; i < FFT_SEGMENT; i++)
				{
					for (int k = 0; k < effLen; k++)
					{
						for (int j = 0; j < TAP_NUM; j++)
						{
							pre_Filter_rx[packet_time[nuser][t] + preamble_len + effLen * i + k][0] += pre_packet[preamble_len + effLen * i + k - j * Unit][0] * h[0][nuser][i][j][0];
							pre_Filter_rx[packet_time[nuser][t] + preamble_len + effLen * i + k][0] -= pre_packet[preamble_len + effLen * i + k - j * Unit][1] * h[0][nuser][i][j][1];
							pre_Filter_rx[packet_time[nuser][t] + preamble_len + effLen * i + k][1] += pre_packet[preamble_len + effLen * i + k - j * Unit][0] * h[0][nuser][i][j][1];
							pre_Filter_rx[packet_time[nuser][t] + preamble_len + effLen * i + k][1] += pre_packet[preamble_len + effLen * i + k - j * Unit][1] * h[0][nuser][i][j][0];
						}
					}
				}

				//----fast fading----
			}
		}
	}

	if (CH_TYPE == 1)
	{
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			if (packet_num[nuser] == 0)
				break;

			for (int i = 0; i < packet_dur * TEST_FRAME_SIZE; i++)
			{
				if (i < TAP_NUM * Unit)
				{
					for (int j = 0; j < i / Unit; j++)
					{
						pre_Filter_rx[i][0] += packet[nuser][i - j * Unit][0] * h[0][nuser][0][j][0];
						pre_Filter_rx[i][0] -= packet[nuser][i - j * Unit][1] * h[0][nuser][0][j][1];
						pre_Filter_rx[i][1] += packet[nuser][i - j * Unit][0] * h[0][nuser][0][j][1];
						pre_Filter_rx[i][1] += packet[nuser][i - j * Unit][1] * h[0][nuser][0][j][0];
					}
				}
				else
				{
					for (int j = 0; j < TAP_NUM; j++)
					{
						pre_Filter_rx[i][0] += packet[nuser][i - j * Unit][0] * h[0][nuser][0][j][0];
						pre_Filter_rx[i][0] -= packet[nuser][i - j * Unit][1] * h[0][nuser][0][j][1];
						pre_Filter_rx[i][1] += packet[nuser][i - j * Unit][0] * h[0][nuser][0][j][1];
						pre_Filter_rx[i][1] += packet[nuser][i - j * Unit][1] * h[0][nuser][0][j][0];
					}
				}
			}
		}
	}
	
	//--- rx Filter
	
	if (PULSE_SHAPE)
	{
		for (int i = 0; i < packet_dur * TEST_FRAME_SIZE; i++)
		{
			for (int j = -trunLen; j < trunLen + UP_RATE; j++)
			{
				rx[i][0] += txFilter[j + trunLen] * pre_Filter_rx[(packet_dur * TEST_FRAME_SIZE + (i - j)) % (packet_dur * TEST_FRAME_SIZE)][0]; // real
				rx[i][1] += txFilter[j + trunLen] * pre_Filter_rx[(packet_dur * TEST_FRAME_SIZE + (i - j)) % (packet_dur * TEST_FRAME_SIZE)][1]; // imaginary
			}
		}
	}
	else
	{
		for (int i = 0; i < packet_dur * TEST_FRAME_SIZE; i++)
		{
			rx[i][0] += pre_Filter_rx[i][0];
			rx[i][1] += pre_Filter_rx[i][1];
		}
	}
		
	
	/*for (int i = 0; i < TEST_FRAME_SIZE * packet_dur; i++)
	{
		cout << rx[i][0] << " ";
	}
	system("pause");
	for (int i = 0; i < TEST_FRAME_SIZE * (FFT_SEGMENT + DIFF_ENC) * effLen; i++)
	{
		cout << rx[i][1] << " ";
	}
	system("pause");*/
	
}
