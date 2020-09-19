#include <iostream>
#include "parameters.h"
#include <cstring>
#include <math.h>
using namespace std;
#pragma warning(disable:4996)

int main()
{
	//---------- check for overflow ----------
	if ((long double)NUM_USER*BLOCK_NUM*BLOCK_LEN < 0)
	{
		cout << "OVERFLOW" << endl;
		system("pause");
		return 0;
	}
	//---------- memory allocation ----------
	int **data = new int*[NUM_USER*NUM_TX];
	double **app = new double*[2 * BLOCK_LEN];
	double **appLlr = new double *[NUM_USER*NUM_TX];
	double **chCoef = new double *[NUM_USER*NUM_TX];
	for (int i = 0; i < NUM_USER*NUM_TX; i++)
	{
		data[i] = new int[BLOCK_LEN];
		appLlr[i] = new double[2 * BLOCK_LEN];
		chCoef[i] = new double[4];
	}
	for (int i = 0; i < 2 * BLOCK_LEN; i++)
		app[i] = new double[NUM_LEVEL];

    double ***supLevel = new double**[NUM_USER*NUM_TX];
    for (int i = 0; i < NUM_USER*NUM_TX; i++){
		supLevel[i] = new double*[NUM_USER * NUM_TX - 1];
		for (int j = 0; j < (NUM_USER * NUM_TX - 1); j++){
				supLevel[i][j] = new double[2]; // real and imaginary
		}
	}


	double** tx = new double* [NUM_USER*NUM_TX];
	double *preTx = nullptr, **rx = nullptr, ****postRx = new double***[NUM_USER], *txFilter = nullptr, *rxFilter = nullptr, *prePilot=nullptr, ** pilot=nullptr;
	
	
	for (int i = 0; i < NUM_USER; i++){
		postRx[i] = new double**[NUM_TX];
		for (int j = 0; j < NUM_TX; j++){
			postRx[i][j] = new double*[BLOCK_LEN];
			for (int k = 0; k < BLOCK_LEN; k++){
			    postRx[i][j][k] = new double[2];
			}
		}
	}
	if (SYNCHRONOUS)
	{
		for (int i = 0; i < NUM_USER*NUM_TX; i++)
		{
			tx[i] = new double[2 * BLOCK_LEN];
		}
		rx = new double *[2 * BLOCK_LEN];
		for (int i = 0; i < 2 * BLOCK_LEN; i++)
		{
			rx[i] = new double[2]; // real and imaginary
		}
	}
	else
	{
		pilot = new double* [NUM_USER*NUM_TX];
		for (int i = 0; i < NUM_USER*NUM_TX; i++)
		{
			tx[i] = new double[BLOCK_LEN*UP_RATE*SPREAD_LEN];
			pilot[i] = new double[BLOCK_LEN*UP_RATE*SPREAD_LEN];
		}
		preTx = new double[BLOCK_LEN*UP_RATE*SPREAD_LEN];
		prePilot = new double[BLOCK_LEN*UP_RATE*SPREAD_LEN];
		rx = new double *[BLOCK_LEN*UP_RATE*SPREAD_LEN];
		postRx = new double*[BLOCK_LEN*UP_RATE*SPREAD_LEN];
		for (int i = 0; i < BLOCK_LEN*UP_RATE*SPREAD_LEN; i++)
		{
			rx[i] = new double[2]; // real and imaginary
			postRx[i] = new double[2]; // real and imaginary
		}
		txFilter = new double[(2 * TRUNCATION + 1)*UP_RATE];
		rxFilter = new double[(2 * TRUNCATION + 1)*UP_RATE];
		SRRCGeneration(rxFilter, 0);
	}
	int *groupSize = nullptr, **group = nullptr;
	double *variation = nullptr, *distList = nullptr, **centroid = nullptr, **estimate = nullptr, **softAssign = nullptr;
	if (CE_SCHEME == 0)
	{
		groupSize = new int[GROUP_SIZE];
		variation = new double[GROUP_SIZE];
		distList = new double[2 * BLOCK_LEN];
		group = new int*[GROUP_SIZE];
		centroid = new double *[GROUP_SIZE];
		for (int i = 0; i < GROUP_SIZE; i++)
		{
			group[i] = new int[2 * BLOCK_LEN];
			centroid[i] = new double[2]; // real and imaginary
		}
		estimate = new double *[NUM_USER*NUM_TX];
		for (int i = 0; i < NUM_USER*NUM_TX; i++)
		{
			estimate[i] = new double[2]; // real and imaginary
		}
		if (EM_GMM)
		{
			softAssign = new double*[2 * BLOCK_LEN];
			for (int i = 0; i < 2 * BLOCK_LEN; i++)
			{
				softAssign[i] = new double[GROUP_SIZE];
			}
		}
	}

	int* known_drift = nullptr;
	//if (COLLISION == 0)
	known_drift = new int[2 * NUM_USER*NUM_TX];
    
	int ***trellis = nullptr;
	double **alpha = nullptr, **beta = nullptr, ***gamma = nullptr;
	if (DIFF_ENC && (DIFF_RX_SCHEME == 1))
	{
		alpha = new double*[sizeof(double)*(BLOCK_LEN + 1)];
		beta = new double*[sizeof(double)*(BLOCK_LEN + 1)];
		for (int i = 0; i < BLOCK_LEN + 1; i++)
		{
			alpha[i] = new double[2];
			beta[i] = new double[2];
		}
		gamma = new double **[BLOCK_LEN];
		for (int i = 0; i < BLOCK_LEN; i++)
		{
			gamma[i] = new double*[2];
			for (int j = 0; j < 2; j++)
			{
				gamma[i][j] = new double[2];
			}
		}
		trellis = new int**[2];
		for (int i = 0; i < 2; i++)
		{
			trellis[i] = new int*[2];
			for (int j = 0; j < 2; j++)
			{
				trellis[i][j] = new int[2];
			}
		}
		//TrellisConstruction(trellis);
	}
	
	double ber[SNR_NUM],bler[SNR_NUM], mse[SNR_NUM] = { 0 }, mse_centroid[SNR_NUM] = { 0 };
	long double itCount[SNR_NUM] = { 0 };
	FILE *result_txt = fopen("result.txt", "w");
	//---------- declaration ----------
	printf("alamouti GDMA-BPSK system\n");
	printf("Number of users: %d\n\n", NUM_USER);
	printf("Number of TXs: %d\n\n", NUM_TX);
	fprintf(result_txt, "GDMA-BPSK system\n");
	fprintf(result_txt, "Number of users: %d\n\n", NUM_USER);
	if (CE_SCHEME == 1)
	{
		printf("Perfect CSI\n\n");
		fprintf(result_txt, "Perfect CSI\n\n");
	}
	else
	{
		if (PROPOSAL == 0)
		{
			printf("Clustering method : original\n");
			fprintf(result_txt, "Clustering method : original\n\n");
		}
		else if (PROPOSAL == 1)
		{
			printf("Clustering method : propose-1\n");
			fprintf(result_txt, "Clustering method : propose-1\n\n");
		}
		else
		{
			printf("Clustering method : propose-2\n");
			fprintf(result_txt, "Clustering method : propose-2\n\n");
		}
		if (INI_METHOD == 0)
		{
			printf("Cluster-based channel estimation: LBG algo.\n");
			fprintf(result_txt, "Cluster-based channel estimation: LBG algo.\n");
		}
		if (INI_METHOD == 1)
		{
			printf("Cluster-based channel estimation: k-means++\n");
			fprintf(result_txt, "Cluster-based channel estimation: k-means++\n");
		}
		if (INI_METHOD == 2)
		{
			printf("Cluster-based channel estimation: modified k-means++\n");
			fprintf(result_txt, "Cluster-based channel estimation: modified k-means++\n");
		}
		if (EM_GMM)
		{
			printf("w/ Gaussian mixture model\n");
			fprintf(result_txt, "w/ Gaussian mixture model\n");
		}
		printf("\n"); fprintf(result_txt, "\n");
	}
	printf("Block length: %d\n\n", BLOCK_LEN);
	fprintf(result_txt, "Block length: %d\n\n", BLOCK_LEN);
	if (DIFF_ENC)
	{
		printf("Differential encoding,\n");
		fprintf(result_txt, "Differential encoding,\n");
	}
	if (DIFF_ENC && (DIFF_RX_SCHEME == 1))
	{
		printf("w/ soft demapping\n\n");
		fprintf(result_txt, "w/ soft demapping\n\n");
	}
	if (DIFF_ENC && (DIFF_RX_SCHEME == 0))
	{
		printf("w/ differential decoding\n\n");
		fprintf(result_txt, "w/ differential decoding\n\n");
	}
	if (SYNCHRONOUS)
	{
		printf("Synchronous system\n\n");
		fprintf(result_txt, "Synchronous system\n\n");
	}
	else
	{
		printf("Asynchronous system\n");
		printf("Drifting range: [%.3f, %.3f]\n\n", -0.5*DRIFT_RANGE, 0.5*DRIFT_RANGE);
		fprintf(result_txt, "Asynchronous system\n");
		fprintf(result_txt, "Drifting range: [%.3f, %.3f]\n\n", -0.5*DRIFT_RANGE, 0.5*DRIFT_RANGE);
	}
	printf("-\n\n"); fprintf(result_txt, "-\n\n");

	//---------- simulation process ----------
	for (int i = 0; i < SNR_NUM; i++)
	{
		double snrdB = SNR_START + (double)i * SNR_STEP;
		double snr = pow(10., snrdB / 10.);
		double variance = 1 / snr;
		double stdDev = sqrt(1 / snr);
		long double error = 0;
		printf("SNR[dB] = %.1f\n", snrdB);
		for (int block = 1; block <= BLOCK_NUM; block++)
		{
			//TODO
			AlamoutiEncoder(data, tx);
			EnergyProfile(chCoef);
			MultipleAccessChannel(sqrt(variance), chCoef, tx, rx, pilot,known_drift);
			if (CE_SCHEME == 0) // cluster-based channel estimation
			{
				Clustering(rx, centroid, group, groupSize, distList, variation, softAssign, variance, estimate, itCount[i],chCoef,known_drift);
				MSEComparison(chCoef, estimate, mse[i],rx,centroid,known_drift); // user specification
				//CentroidMSEComparison(chCoef, estimate, mse_centroid[i], rx, centroid);
			}
			SignalCombiner(chCoef, rx, postRx);
			SuperLevelSpecification(chCoef, supLevel);
			MLDT(pow(stdDev, 2), chCoef, supLevel, postRx, app, appLlr);





			Detector(data, appLlr, trellis, alpha, beta, gamma, error,known_drift);
			ber[i] = error / ((long double)NUM_USER * NUM_TX * block * (BLOCK_LEN - DIFF_ENC));
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