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

void EnergyProfile(double ****h, double ****H)
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
			for (int j = 0; j < FFT_POINT; j++)
			{

				if (j < TAP_NUM)
				{
					double power = 0.5*exp(-j) / nFactor;
					double phi = uniform(generator) * 2. * M_PI - M_PI;
					h[i][0][j][0] = normal(generator);
					h[i][0][j][1] = normal(generator);
					h[i][0][j][0] = sqrt(power) * (pow(K + 1, -0.5) * h[i][0][j][0] + pow(K, 0.5) / pow(K + 1, 0.5) * cos(phi));
					h[i][0][j][1] = sqrt(power) * (pow(K + 1, -0.5) * h[i][0][j][1] + pow(K, 0.5) / pow(K + 1, 0.5) * sin(phi));
					//cout << h[i][0][j][0] << " " << h[i][0][j][1] << endl;
				}
				else h[i][0][j][0] = h[i][0][j][1] = 0;
				if (i % 2 == 1){
					h[i][0][j][0] = h[i-1][0][j][1];
					h[i][0][j][1] = (-1.0)*h[i-1][0][j][0];
				}
				for (int k = 1; k < FFT_SEGMENT + DIFF_ENC; k++) // block fading
				{
					h[i][k][j][0] = h[i][0][j][0];
					h[i][k][j][1] = h[i][0][j][1];
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
					double power = 0.5*exp(-j) / nFactor;
					double theta = uniform(generator) * 2. * M_PI - M_PI;
					for (int k = 0; k < FFT_SEGMENT + DIFF_ENC; k++)
					{
						h[i][k][j][0] = h[i][k][j][1] = 0;
					}
					// Derive c_i and c_q
					
					for (int n = 1; n <= M; n++)
					{
						double phi_1 = uniform(generator) * 2. * M_PI - M_PI;
						double phi_2 = uniform(generator) * 2. * M_PI - M_PI;
						double alpha = (2.*M_PI*n - M_PI + theta) / (4.*M);
						for (int k = 0; k < FFT_SEGMENT + DIFF_ENC; k++)
						{
							double t = k*FFT_POINT*period;
							h[i][k][j][0] += cos(2.*M_PI*fd*t*cos(alpha) + phi_1);
							h[i][k][j][1] += cos(2.*M_PI*fd*t*sin(alpha) + phi_2);
						}
					}
					// Derive different channel tap
					double phi = uniform(generator) * 2. * M_PI - M_PI;
					double thi = uniform(generator) * 2. * M_PI - M_PI;
					for (int k = 0; k < FFT_SEGMENT + DIFF_ENC; k++)
					{
						double t = k * FFT_POINT * period;
						h[i][k][j][0] = sqrt(power) * (pow(K + 1, -0.5) * h[i][k][j][0] * sqrt(2. / M) + pow(K, 0.5) / pow(K + 1, 0.5) * cos(2. * M_PI * fd * t * cos(thi) + phi));
						h[i][k][j][1] = sqrt(power) * (pow(K + 1, -0.5) * h[i][k][j][1] * sqrt(2. / M) + pow(K, 0.5) / pow(K + 1, 0.5) * sin(2. * M_PI * fd * t * cos(thi) + phi));
					}
				}
				else
				{
					// Derive different channel tap
					for (int k = 0; k < FFT_SEGMENT + DIFF_ENC; k++)
					{
						h[i][k][j][0] = h[i][k][j][1] = 0;
					}
				}
			}
		}
	}

	//printchannel(h);
	//---------- channel coefficients in the frequency domain ----------
	double real[FFT_POINT], imag[FFT_POINT];
	for (int i = 0; i <NUM_USER; i++)
	{
		for (int j = 0; j < FFT_SEGMENT + DIFF_ENC; j++)
		{
			for (int k = 0; k < FFT_POINT; k++)
			{
				real[k] = h[i][j][k][0];
				imag[k] = h[i][j][k][1];
			}
			FFT(real, imag, FFT_POINT, FFT_LAYER, H[i][j][0], H[i][j][1], 0, 0);
		}
	}
}

void MultipleAccessChannel(double stdDev, double ****h, double ****tx, double ***rx, double *drift, double ****H)
{
	int effLen; //-------------�@��block�T��������
	int effTap; //-------------�@��symbol�T��������

	// effLen & effTap setting
	if (SYNCHRONOUS)
	{
		effLen = CP_TYPE ? FFT_POINT + CP_LEN : FFT_POINT + CP_LEN + CS_LEN;
		effTap = 1;
	}
	else
	{
		if (CP_TYPE)
		{
			if (OVER_SAMPLE)
			{
				if (PULSE_SHAPE)
				{
					effLen = (FFT_POINT + CP_LEN) * UP_RATE * OVER_SAMPLE_RATE;
					effTap = UP_RATE * OVER_SAMPLE_RATE;
				}
				else
				{
					effLen = (FFT_POINT + CP_LEN) * OVER_SAMPLE_RATE;
					effTap = OVER_SAMPLE_RATE;
				}
			}
			else
			{
				if (PULSE_SHAPE)
				{
					effLen = (FFT_POINT + CP_LEN) * UP_RATE;
					effTap = UP_RATE;
				}
				else
				{
					cout << "~over sample ~pulse shaping ~non syn";
					system("pause");
				}
			}
		}
		else
		{
			if (OVER_SAMPLE)
			{
				if (PULSE_SHAPE)
				{
					effLen = (FFT_POINT + CP_LEN + CS_LEN) * UP_RATE * OVER_SAMPLE_RATE;
					effTap = UP_RATE * OVER_SAMPLE_RATE;
				}
				else
				{
					effLen = (FFT_POINT + CP_LEN + CS_LEN) * OVER_SAMPLE_RATE;
					effTap = OVER_SAMPLE_RATE;
				}
			}
			else
			{
				if (PULSE_SHAPE)
				{
					effLen = (FFT_POINT + CP_LEN + CS_LEN) * UP_RATE;
					effTap = UP_RATE;
				}
				else
				{
					cout << "~over sample ~pulse shaping ~non syn";
					system("pause");
				}
			}
		}
		///////////////////////////////////////

		for (int i = 0; i < NUM_USER; i++)
		{
			drift[i] = uniform_drift(generator);
		}

		//drift[1] -= 4;
		//drift[0] = 0;
		//drift[1] = 0;
		//drift[1] = -5;

		/*cout << endl;
		for (int i = 0; i < effLen; i++)
		{
			cout << tx[0][0][i][0] << " ";
		}*/

		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{

			double drift_ = drift[nuser] * (double)effTap;
			double a = drift_ > 0 ? drift_ - int(drift_) : int(drift_) - drift_;
			double b = drift_ > 0 ? 1 - drift_ + int(drift_) : 1 + drift_ - int(drift_);
			//cout << drift[nuser]<<" "<<drift_<<" "<<int(drift_)<<" "<<a << " " << b << endl;
			for (int i = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
			{
				//----------------���tdrift��
				if (drift_ > 0) // -------------------forword delay
				{
					for (int j = 0; j < effLen; j++)
					{
						if (j < effLen - int(drift_) - 1)
						{
							tx[nuser][i][j][0] = b * tx[nuser][i][j + int(drift_)][0] + a * tx[nuser][i][j + int(drift_) + 1][0];
							tx[nuser][i][j][1] = b * tx[nuser][i][j + int(drift_)][1] + a * tx[nuser][i][j + int(drift_) + 1][1];
						}
						else //-----------------IBI
						{
							if (i != FFT_SEGMENT + DIFF_ENC - 1)
							{
								tx[nuser][i][j][0] = b * tx[nuser][i + 1][(j - (effLen - int(drift_) - 1)) % effLen][0] + a * tx[nuser][i + 1][(j - (effLen - int(drift_) - 1) + 1) % effLen][0];
								tx[nuser][i][j][1] = b * tx[nuser][i + 1][(j - (effLen - int(drift_) - 1)) % effLen][1] + a * tx[nuser][i + 1][(j - (effLen - int(drift_) - 1) + 1) % effLen][1];
							}
							else
							{
								tx[nuser][i][j][0] = 0;
								tx[nuser][i][j][1] = 0;
							}
						}
					}
				}
				else // -------------------backword delay
				{
					for (int j = effLen - 1; j >= 0; j--)
					{
						if (j > -int(drift_))
						{
							tx[nuser][i][j][0] = a * tx[nuser][i][j + int(drift_)][0] + b * tx[nuser][i][j + int(drift_) - 1][0];
							tx[nuser][i][j][1] = a * tx[nuser][i][j + int(drift_)][1] + b * tx[nuser][i][j + int(drift_) - 1][1];
						}
						else  //---------------- IBI
						{
							if (i != 0)
							{
								tx[nuser][i][j][0] = a * tx[nuser][i - 1][(j + int(drift_) + effLen) % effLen][0] + b * tx[nuser][i][(j + int(drift_) - 1 + effLen) % effLen][0];
								tx[nuser][i][j][1] = a * tx[nuser][i - 1][(j + int(drift_) + effLen) % effLen][1] + b * tx[nuser][i][(j + int(drift_) - 1 + effLen) % effLen][1];
							}
							else
							{
								tx[nuser][i][j][0] = 0;
								tx[nuser][i][j][1] = 0;
							}
						}
					}
				}
				/*	cout << endl;
					for (int j = 0; j < effLen; j++)
					{
						cout << tx[nuser][i][j][0] << " ";
					}
					system("pause");*/
			}
		}

		//--- channel coefficient shift
		/*for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			for (int j = 0; j < 2; j++)
			{
				for (int i = 0; i < FFT_POINT; i++)
				{
					cout << H[nuser][0][j][i] << " ";
				}
				cout << endl;
			}
			cout << endl;
		}*/

		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			for (int i = 0; i < FFT_POINT; i++)
			{
				for (int j = 0; j < FFT_SEGMENT + DIFF_ENC; j++)
				{
					vector<double> Phasor(2);

					Phasor[0] = sqrt(pow(H[nuser][j][0][i], 2) + pow(H[nuser][j][1][i], 2)); // amplitude
					Phasor[1] = (H[nuser][j][1][i] > 0) ? (acos(H[nuser][j][0][i] / Phasor[0])) : (2. * M_PI - acos(H[nuser][j][0][i] / Phasor[0])); // phase
					Phasor[1] += 2 * M_PI * i * drift[nuser] / FFT_POINT;
					H[nuser][j][0][i] = Phasor[0] * cos(Phasor[1]);
					H[nuser][j][1][i] = Phasor[0] * sin(Phasor[1]);
				}
			}
		}

		/*for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			for (int j = 0; j < 2; j++)
			{
				for (int i = 0; i < FFT_POINT; i++)
				{
					cout << H[nuser][0][j][i] << " ";
				}
				cout << endl;
			}
			cout << endl;
		}
		system("pause");*/
	}

	// MAC
	for (int i = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
	{
		for (int j = 0; j < effLen; j++)
		{
			rx[i][j][0] = stdDev * normal(generator); // real
			rx[i][j][1] = stdDev * normal(generator); // imaginary
			//rx[i][j][0] = 0; // real
			//rx[i][j][1] = 0; // imaginary
		}
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			for (int j = 0; j < effLen; j++)
			{
				for (int k = 0; k < TAP_NUM; k++)
				{
					rx[i][j][0] += tx[nuser][i][(effLen + (j - k * effTap) % effLen) % effLen][0] * h[nuser][i][k][0];
					rx[i][j][0] -= tx[nuser][i][(effLen + (j - k * effTap) % effLen) % effLen][1] * h[nuser][i][k][1];
					rx[i][j][1] += tx[nuser][i][(effLen + (j - k * effTap) % effLen) % effLen][1] * h[nuser][i][k][0];
					rx[i][j][1] += tx[nuser][i][(effLen + (j - k * effTap) % effLen) % effLen][0] * h[nuser][i][k][1];
				}
			}
		}
	}
}
