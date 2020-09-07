#include <iostream>
#include "parameters.h"
#include <cstring>
#include <math.h>
using namespace std;
#pragma warning(disable:4996)

int main()
{
	//---------- check for overflow ----------
	if ((long double)NUM_USER*NUM_TX*BLOCK_NUM < 0)
	{
		cout << "OVERFLOW" << endl;
		system("pause");
		return 0;
	}
	//---------- memory allocation ----------
	int **data = new int*[NUM_USER];
	double **appLlr = new double *[NUM_USER];
	double ***tx = new double**[NUM_USER];
	double ***chCoef = new double **[NUM_USER];
	double ***postRx = new double**[NUM_USER];
	double ***app = new double**[NUM_USER];
	double ****supLevel = new double***[NUM_USER];
	for (int i = 0; i < NUM_USER; i++)
	{
		data[i] = new int[L];
		appLlr[i] = new double[NUM_TX];
		tx[i] = new double*[L];
		chCoef[i] = new double*[NUM_TX];
		postRx[i] = new double*[L];
		app[i] = new double*[NUM_TX];
		supLevel[i] = new double**[NUM_TX];
		for (int j = 0; j < NUM_TX; j++)
		{
			chCoef[i][j] = new double[2]; // real and imaginary
			app[i][j] = new double[NUM_LEVEL];
			supLevel[i][j] = new double*[NUM_USER * NUM_TX - 1];
			for (int k = 0; k < (NUM_USER * NUM_TX - 1); k++)
			{
				supLevel[i][j][k] = new double[2]; // real and imaginary
			}
		}
		for (int j = 0; j < L; j++)
		{
			tx[i][j] = new double[NUM_TX];
			postRx[i][j] = new double[2]; // real and imaginary
		}
	}
	double **rx = new double*[L];
	for (int i = 0; i < L; i++)
	{
		rx[i] = new double[2]; // real and imaginary
	}
	double ber[SNR_NUM];
	int** STBC_matrix = new int*[NUM_TX];
	for (int i = 0; i < NUM_TX; i++)
	{
		STBC_matrix[i] = new int[L];
	}
	int** SC_matrix = new int* [L];
	for (int i = 0; i < L; i++)
	{
		SC_matrix[i] = new int[NUM_TX];
	}
	FILE *result_txt = fopen("result.txt", "w");
	//---------- declaration ----------
	printf("STC-GDMA-BPSK system\n");
	printf("Number of users: %d\n\n", NUM_USER);
	fprintf(result_txt, "STC-GDMA-BPSK system\n");
	fprintf(result_txt, "Number of users: %d\n\n", NUM_USER);
	printf("-\n\n"); fprintf(result_txt, "-\n\n");
	//---------- simulation process ----------
	for (int i = 0; i < SNR_NUM; i++)
	{
		double snrdB = SNR_START + (double)i * SNR_STEP;
		double snr = pow(10., snrdB / 10.);
		double stdDev = sqrt(1 / snr);
		long double error = 0;
		printf("SNR[dB] = %.1f\n", snrdB);
		for (int block = 1; block <= BLOCK_NUM; block++)
		{
			AlamoutiEncoder(data, tx, STBC_matrix);
			EnergyProfile(chCoef);
			MultipleAccessChannel(stdDev, chCoef, tx, rx);
			SignalCombiner(chCoef, rx, postRx, SC_matrix);
			SuperLevelSpecification(chCoef, supLevel, SC_matrix);
			MLDT(pow(stdDev, 2), chCoef, supLevel, postRx, app, appLlr);
			Detector(data, appLlr, error);
			ber[i] = error / ((long double)NUM_USER * NUM_TX * block);
			printf("Block# = %d, BER = %e\r", block, ber[i]);
		}
		fprintf(result_txt, "SNR[dB] = %.1f, BER = %e\n", snrdB, ber[i]);
		printf("\n");
	}
	fprintf(result_txt, "\n\n");
	for (int i = 0; i < SNR_NUM; i++)
	{
		fprintf(result_txt, "%e, ", ber[i]);
	}
	fprintf(result_txt, "\n");
	fclose(result_txt);
	system("pause");
	return 0;
}