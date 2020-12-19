#include <iostream>
#include "parameters.h"
#include <cstring>
#include <math.h>
using namespace std;
#pragma warning(disable:4996)

//--------- scma matrix-------------//can be generalize
int		 scma_matrix[SCMA_SOURCE][SCMA_SOURCE * NUM_USER / SCMA_USER_SOURCE] = { {0,1,1,0,1,0}, {1,0,1,0,0,1}, {0,1,0,1,0,1}, {1,0,0,1,1,0} };
int		 coef_idx[SCMA_SOURCE][SCMA_SOURCE * NUM_USER / SCMA_USER_SOURCE] = { {-1,0,0,-1,0,-1}, {0,-1,1,-1,-1,0}, {-1,1,-1,0,-1,1}, {1,-1,-1,1,1,-1} };

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
    int ***data = new int**[SCMA_SOURCE];
	double ***appLlr = new double **[SCMA_SOURCE];
	double ****tx = new double***[SCMA_SOURCE];
	double ****chCoef = new double ***[SCMA_SOURCE];
	double ****postRx = new double***[SCMA_SOURCE];
	double ****app = new double***[SCMA_SOURCE];
	double *****supLevel = new double****[SCMA_SOURCE];
	for (int q = 0; q < SCMA_SOURCE;q++){
        data[q] = new int*[NUM_USER];
	    appLlr[q] = new double *[NUM_USER];
	    tx[q] = new double**[NUM_USER];
	    chCoef[q] = new double **[NUM_USER];
	    postRx[q] = new double**[NUM_USER];
	    app[q] = new double**[NUM_USER];
	    supLevel[q] = new double***[NUM_USER];
		for (int i = 0; i < NUM_USER; i++)
		{
			data[q][i] = new int[NUM_TX];
			appLlr[q][i] = new double[NUM_TX];
			tx[q][i] = new double*[NUM_TX];
			chCoef[q][i] = new double*[NUM_TX];
			postRx[q][i] = new double*[NUM_TX];
			app[q][i] = new double*[NUM_TX];
			supLevel[q][i] = new double**[NUM_TX];
			for (int j = 0; j < NUM_TX; j++)
			{
				tx[q][i][j] = new double[NUM_TX];
				chCoef[q][i][j] = new double[2]; // real and imaginary
				postRx[q][i][j] = new double[2]; // real and imaginary
				app[q][i][j] = new double[NUM_LEVEL];
				supLevel[q][i][j] = new double*[NUM_USER * NUM_TX - 1];
				for (int k = 0; k < (NUM_USER * NUM_TX - 1); k++)
				{
					supLevel[q][i][j][k] = new double[2]; // real and imaginary
				}
			}
		}
	}
	double*** rx = new double** [SCMA_SOURCE];
	for (int q = 0; q < SCMA_SOURCE;q++){
        rx[q] = new double*[NUM_TX];
		for (int i = 0; i < NUM_TX; i++)
	    {
	    	rx[q][i] = new double[2]; // real and imaginary
	    }
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
			for (int q = 0; q < SCMA_SOURCE;q++){
			    AlamoutiEncoder(data[q], tx[q]);//TODO:send same message
			    EnergyProfile(chCoef[q]);
			    MultipleAccessChannel(stdDev, chCoef[q], tx[q], rx[q]);
			    SignalCombiner(chCoef[q], rx[q], postRx[q]);
			    SuperLevelSpecification(chCoef[q], supLevel[q]);
				MLDT(pow(stdDev, 2), chCoef[q], supLevel[q], postRx[q], app[q], appLlr[q]);

			}
//////TODO!!! and recalculate stdDev
			for (int q = 0; q < SCMA_SOURCE;q++){

			    Detector(data[q], appLlr[q], error);
			}
			ber[i] = error / ((long double)NUM_USER * NUM_TX * block * SCMA_SOURCE);
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
	return 0;
}