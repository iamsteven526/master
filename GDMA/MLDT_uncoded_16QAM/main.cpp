#include <iostream>
#include "parameters.h"
#include <cmath>
#include <vector>
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
    //16QAM
    vector<vector<double>> map_16QAM(16, vector<double> (2, 0));

    map_16QAM[0][0] = -0.94868;
	map_16QAM[0][1] = 0.94868;

	map_16QAM[1][0] = -0.94868;
	map_16QAM[1][1] = -0.94868;

	map_16QAM[2][0] = 0.94868;
	map_16QAM[2][1] = 0.94868;

    map_16QAM[3][0] = 0.94868;
	map_16QAM[3][1] = -0.94868;

	map_16QAM[4][0] = -0.94868;
	map_16QAM[4][1] = 0.31622;

	map_16QAM[5][0] = -0.31622;
	map_16QAM[5][1] = -0.94868;

    map_16QAM[6][0] = 0.31622;
	map_16QAM[6][1] = 0.94868;

	map_16QAM[7][0] = 0.94868;
	map_16QAM[7][1] = -0.31622;

	map_16QAM[8][0] = -0.31622;
	map_16QAM[8][1] = 0.94868;

    map_16QAM[9][0] = -0.94868;
	map_16QAM[9][1] = -0.31622;

	map_16QAM[10][0] = 0.94868;
	map_16QAM[10][1] = 0.31622;

	map_16QAM[11][0] = 0.31622;
	map_16QAM[11][1] = -0.94868;

    map_16QAM[12][0] = -0.31622;
	map_16QAM[12][1] = 0.31622;

	map_16QAM[13][0] = -0.31622;
	map_16QAM[13][1] = -0.31622;

	map_16QAM[14][0] = 0.31622;
	map_16QAM[14][1] = 0.31622;

	map_16QAM[15][0] = 0.31622;
	map_16QAM[15][1] = -0.31622;





	//---------- memory allocation ----------
	int **data = new int*[NUM_USER];
	double **app = new double*[2 * BLOCK_LEN];
	double **appLlr = new double *[NUM_USER];
	double **chCoef = new double *[NUM_USER];
	for (int i = 0; i < NUM_USER; i++)
	{
		data[i] = new int[BLOCK_LEN];
		appLlr[i] = new double[2 * BLOCK_LEN];
		chCoef[i] = new double[4];
	}
	for (int i = 0; i < 2 * BLOCK_LEN; i++)
		app[i] = new double[NUM_LEVEL];

	double*** tx = new double** [NUM_USER];
	double *preTx = nullptr, **rx = nullptr, **postRx = nullptr, *txFilter = nullptr, *rxFilter = nullptr, *prePilot=nullptr, ** pilot=nullptr;
	if (SYNCHRONOUS)
	{
		for (int i = 0; i < NUM_USER; i++)
		{
			tx[i] = new double *[2 * BLOCK_LEN / MOD_LEVEL];
				for (int j = 0; j < 2 * BLOCK_LEN / MOD_LEVEL; j++)
				{
					tx[i][j] = new double[2]; // real and imaginary
				}
		}
		rx = new double *[2 * BLOCK_LEN / MOD_LEVEL];
		for (int i = 0; i < 2 * BLOCK_LEN / MOD_LEVEL; i++)
		{
			rx[i] = new double[2]; // real and imaginary
		}
	}
	/*
	else
	{
		pilot = new double* [NUM_USER];
		for (int i = 0; i < NUM_USER; i++)
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
	*/
	int *groupSize = nullptr, **group = nullptr;
	double *variation = nullptr, *distList = nullptr, **centroid = nullptr, **estimate = nullptr, **softAssign = nullptr;
	if (CE_SCHEME == 0)
	{
		groupSize = new int[GROUP_SIZE];
		variation = new double[GROUP_SIZE];
		distList = new double[2 * BLOCK_LEN / MOD_LEVEL];
		group = new int*[GROUP_SIZE];
		centroid = new double *[GROUP_SIZE];
		for (int i = 0; i < GROUP_SIZE; i++)
		{
			group[i] = new int[2 * BLOCK_LEN / MOD_LEVEL];
			centroid[i] = new double[2]; // real and imaginary
		}
		estimate = new double *[NUM_USER];
		for (int i = 0; i < NUM_USER; i++)
		{
			estimate[i] = new double[2]; // real and imaginary
		}
		if (EM_GMM)
		{
			softAssign = new double*[2 * BLOCK_LEN / MOD_LEVEL];
			for (int i = 0; i < 2 * BLOCK_LEN / MOD_LEVEL; i++)
			{
				softAssign[i] = new double[GROUP_SIZE];
			}
		}
	}

	int* known_drift = nullptr;
	//if (COLLISION == 0)
	known_drift = new int[2 * NUM_USER];

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
		TrellisConstruction(trellis);
	}
	double ber[SNR_NUM],bler[SNR_NUM], mse[SNR_NUM] = { 0 }, mse_centroid[SNR_NUM] = { 0 };
	long double itCount[SNR_NUM] = { 0 };
	FILE *result_txt = fopen("result.txt", "w");
	//---------- declaration ----------
	printf("GDMA-BPSK system\n");
	printf("Number of users: %d\n\n", NUM_USER);
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
		double variance = 0.5 / (4*snr);
		long double error = 0;
		printf("SNR[dB] = %.1f\n", snrdB);
		int block = 1;
		for (block; block <= BLOCK_NUM; block++)
		{
			
			//cout << "meow" << endl;
			Transmitter(data, preTx, tx, txFilter,rxFilter,pilot,prePilot, known_drift);
			//Transmitter(data, preTx, tx, txFilter,pilot);
			EnergyProfile(chCoef);
			MultipleAccessChannel(sqrt(variance), chCoef, tx, rx, pilot,known_drift);
			if (!SYNCHRONOUS) ReceivingFilter(rx, postRx, rxFilter, txFilter, pilot,chCoef);
			//system("pause");
			if (CE_SCHEME == 0) // cluster-based channel estimation
			{
				Clustering(rx, centroid, group, groupSize, distList, variation, softAssign, variance, estimate, itCount[i],chCoef,known_drift);
				MSEComparison(chCoef, estimate, mse[i],rx,centroid,known_drift); // user specification
				//CentroidMSEComparison(chCoef, estimate, mse_centroid[i], rx, centroid);
			}
			MLDT(variance, chCoef, rx, app, appLlr, estimate,known_drift, map_16QAM);
			int BLERa = error;
			Detector(data, appLlr, trellis, alpha, beta, gamma, error,known_drift);
			int BLERb = error;
            if(BLERa != BLERb){
				bler[i]++;
			}
			ber[i] = error / ((long double)NUM_USER * block * (BLOCK_LEN - DIFF_ENC)); // excluding reference symbol
			//cout << itCount[i] << " ";
			if (CE_SCHEME == 1) printf("Block# = %d, BER = %e\r", block, ber[i]); 
			else printf("Block# = %d, BER = %e,BLER = %e, MSE = %e, MSE_centroid = %e,iteration = %e\r", block, ber[i],(bler[i]/block), mse[i] / (double)(NUM_USER * block), mse_centroid[i] / (double)(NUM_LEVEL * block), itCount[i]/(double)(block));
			
			if (error > 10000 && block>500000)
				break;
		}
		mse[i] /= (double)(NUM_USER* block);
		itCount[i] /= (double)(NUM_USER * block);
		if (CE_SCHEME == 1) fprintf(result_txt, "SNR[dB] = %.1f, BER = %e\n", snrdB, ber[i]); 
		else fprintf(result_txt, "SNR[dB] = %.1f, BER = %e, MSE = %e, MSE_centroid = %e\n", snrdB, ber[i], mse[i], mse_centroid[i]);
		printf("\n\n");
	}
	fprintf(result_txt, "\n\n");
	for (int i = 0; i < SNR_NUM; i++)
	{
		fprintf(result_txt, "%e, ", ber[i]);
	}
	fprintf(result_txt, "\n");
	if (CE_SCHEME == 0)
	{
		for (int i = 0; i < SNR_NUM; i++)
		{
			fprintf(result_txt, "%e, ", mse[i]);
		}
		fprintf(result_txt, "\n");
		for (int i = 0; i < SNR_NUM; i++)
		{
			fprintf(result_txt, "%f, ", itCount[i]);
		}
		fprintf(result_txt, "\n");
	}
	fclose(result_txt);
	//system("pause");
	return 0;
}
