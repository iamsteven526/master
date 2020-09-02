#include <iostream>
#include "parameters.h"
using namespace std;
#pragma warning (disable : 4996)

int main()
{
	//---------- check for overflow ----------
	if ((long double)NUM_USER * BLOCK_NUM * BLOCK_LEN * Qm < 0)
	{
		cout << "OVERFLOW" << endl;
		system("pause");
		return 0;
	}
	//---------- memory allocation ----------
	int** data = new int* [NUM_USER];
	double* app = new double[NUM_LEVEL];
	double** tx = new double* [NUM_USER];
	double** appLlr = new double* [NUM_USER];
	double** chCoef = new double* [NUM_USER];
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		data[nuser] = new int[BLOCK_LEN * Qm];
		tx[nuser] = new double[BLOCK_LEN * Qm];
		appLlr[nuser] = new double[BLOCK_LEN * Qm];
		chCoef[nuser] = new double[5];
	}
	double** rx = new double* [BLOCK_LEN];
	for (int i = 0; i < BLOCK_LEN; i++)
	{
		rx[i] = new double[2]; // real and imaginary
	}

	int* groupSize = nullptr, ** group = nullptr;
	double* variation = nullptr, * distList = nullptr, ** centroid = nullptr, ** estimate = nullptr, ** softAssign = nullptr;
	if (CE_METHOD == 0)
	{
		groupSize = new int[GROUP_SIZE];
		variation = new double[GROUP_SIZE];
		distList = new double[2*Cluster_Len];
		group = new int* [GROUP_SIZE];
		centroid = new double* [GROUP_SIZE];
		for (int i = 0; i < GROUP_SIZE; i++)
		{
			group[i] = new int[BLOCK_LEN];
			centroid[i] = new double[2]; // real and imaginary
		}
		estimate = new double* [NUM_USER];
		for (int i = 0; i < NUM_USER; i++)
		{
			estimate[i] = new double[2]; // real and imaginary
		}
		if (EM_GMM)
		{
			softAssign = new double* [BLOCK_LEN];
			for (int i = 0; i < BLOCK_LEN; i++)
			{
				softAssign[i] = new double[GROUP_SIZE];
			}
		}
	}
	FILE* result_txt = fopen("result.txt", "w");
	//---------- declaration ----------
	printf("Uncoded GDMA-QPSK system\n");
	printf("Number of users: %d\n\n", NUM_USER);
	printf("Block length: %d\n\n", BLOCK_LEN);
	printf("Cluster length: %d\n\n", Cluster_Len);
	fprintf(result_txt, "GDMA-QPSK system\n");
	fprintf(result_txt, "Number of users: %d\n\n", NUM_USER);
	fprintf(result_txt, "Block length: %d\n\n", BLOCK_LEN);
	fprintf(result_txt, "Cluster length: %d\n\n",Cluster_Len);
	if (CE_METHOD == 1)
	{
		printf("Perfect CSI\n\n");
		fprintf(result_txt, "Perfect CSI\n\n");
	}
	else
	{
		if (INI_METHOD == 0)
		{
			printf("Cluster-based channel estimation: K-means++\n");
			fprintf(result_txt, "Cluster-based channel estimation: K-means++\n");
		}
		if (INI_METHOD == 1)
		{
			printf("Cluster-based channel estimation: LBG algo.\n");
			fprintf(result_txt, "Cluster-based channel estimation: LBG algo.\n");
		}
		if (INI_METHOD == 2)
		{
			printf("Cluster-based channel estimation: K-means\n");
			fprintf(result_txt, "Cluster-based channel estimation: K-means\n");
		}
		if (INI_METHOD == 3)
		{
			printf("Cluster-based channel estimation: K-means++ added threshold\n");
			fprintf(result_txt, "Cluster-based channel estimation: K-means++\n");
		}
		if (INI_METHOD == 4)
		{
			printf("Cluster-based channel estimation: K-means++ far-point\n");
			fprintf(result_txt, "Cluster-based channel estimation: K-means++\n");
		}
		if (EM_GMM)
		{
			printf("w/Gaussian mixture model\n");
			fprintf(result_txt, "w/Gaussian mixture model\n");
		}
		printf("\n"); fprintf(result_txt, "\n");
	}
	if (DIFF_ENC)
	{
		printf("Differential encoding,\n");
		printf("w/ differential decoding\n\n");
		fprintf(result_txt, "Differential encoding,\n");
		fprintf(result_txt, "w/ differential decoding\n\n");
	}
	double ber[SNR_NUM], mse[SNR_NUM] = { 0 }, mse_centroid[SNR_NUM] = { 0 };
	long double itCount[SNR_NUM] = { 0 };
	printf("-\n\n"); fprintf(result_txt, "-\n\n");
	//---------- simulation process ----------
	for (int i = 0; i < SNR_NUM; i++)
	{
		double snrdB = SNR_START + (double)i * SNR_STEP;
		double snr = pow(10., snrdB / 10.);
		double variance = 0.5 / snr;
		double stdDev = sqrt(0.25 / snr);
		long double error = 0;
		printf("SNR[dB] = %.1f\n", snrdB);
		for (int block = 1; block <= BLOCK_NUM; block++)
		{
			Transmitter(data, tx);
			EnergyProfile(chCoef);
			MultipleAccessChannel(stdDev, chCoef, tx, rx);
			if (CE_METHOD == 0) // cluster-based channel estimation
			{
				KmeansClustering(rx, centroid, group, groupSize, distList, variation, softAssign, variance, estimate, itCount[i], chCoef);
				MSEComparison(chCoef, estimate, mse[i], rx, centroid); // user specification
				CentroidMSEComparison(chCoef, estimate, mse_centroid[i], rx, centroid);
			}
			MLDT(pow(stdDev, 2), chCoef, rx, app, appLlr,estimate);
			Detecter(data, appLlr, error);
			ber[i] = error / ((long double)NUM_USER * block * (BLOCK_LEN - DIFF_ENC) * Qm); // excluding reference symbol
			if (CE_METHOD == 1) printf("Block# = %d, BER = %e\r", block, ber[i]);
			else printf("Block# = %d, BER = %e, MSE = %e, MSE_centroid = %e\r", block, ber[i], mse[i] / (double)(NUM_USER * block), mse_centroid[i] / (double)(NUM_LEVEL * block));
		}
		mse[i] /= (double)(NUM_USER * BLOCK_NUM);
		mse_centroid[i] /= (double)(NUM_LEVEL * BLOCK_NUM);
		itCount[i] /= (double)BLOCK_NUM;
		if (CE_METHOD == 1) fprintf(result_txt, "SNR[dB] = %.1f, BER = %e\n", snrdB, ber[i]);
		else fprintf(result_txt, "SNR[dB] = %.1f, BER = %e, MSE = %e, MSE_centroid = %e\n", snrdB, ber[i], mse[i], mse_centroid[i]);
		printf("\n\n");
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

