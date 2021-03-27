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
		data[i] = new int[NUM_TX];
		appLlr[i] = new double[NUM_TX];
		tx[i] = new double*[NUM_TX];
		chCoef[i] = new double*[NUM_TX];
		postRx[i] = new double*[NUM_TX];
		app[i] = new double*[NUM_TX];
		supLevel[i] = new double**[NUM_TX];
		for (int j = 0; j < NUM_TX; j++)
		{
			tx[i][j] = new double[NUM_TX];
			chCoef[i][j] = new double[2]; // real and imaginary
			postRx[i][j] = new double[2]; // real and imaginary
			app[i][j] = new double[NUM_LEVEL];
			supLevel[i][j] = new double*[NUM_USER * NUM_TX - 1];
			for (int k = 0; k < (NUM_USER * NUM_TX - 1); k++)
			{
				supLevel[i][j][k] = new double[2]; // real and imaginary
			}
		}
	}
	double **rx = new double*[NUM_TX];
	for (int i = 0; i < NUM_TX; i++)
	{
		rx[i] = new double[2]; // real and imaginary
	}
	double ber[SNR_NUM];
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
			AlamoutiEncoder(data, tx);
			EnergyProfile(chCoef);
			MultipleAccessChannel(stdDev, chCoef, tx, rx);
			SignalCombiner(chCoef, rx, postRx);
			SuperLevelSpecification(chCoef, supLevel);
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