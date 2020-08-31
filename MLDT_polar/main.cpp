#define _USE_MATH_DEFINES
#include <iostream>
#include "parameters.h"

using namespace std;

#pragma warning(disable:4996)

int main()
{

	//---------- check for overflow ----------
	if ((long double)BLOCK_NUM*DATA_LEN*NUM_USER < 0)
	{
		cout << "OVERFLOW" << endl;
		system("pause");
		return 0;
	}
	//---------- memory allocation ----------
	vector<vector<uint8_t>> data;
	vector<vector<uint8_t>> codeword;
	vector<vector<double>> appLlr;
	data.resize(NUM_USER);
	codeword.resize(NUM_USER);
	appLlr.resize(NUM_USER);
	double **tx = new double *[NUM_USER];
	double ***chCoef = new double **[NUM_USER];
	double **app = new double *[CODE_LEN];
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		data.at(nuser).resize(DATA_LEN);
		appLlr[nuser].resize(CODE_LEN);
		chCoef[nuser] = new double*[4];
		for (int i = 0; i < 4; i++)
		{
			chCoef[nuser][i] = new double[CODE_LEN/FADING_SIZE+1];
		}
	}
	for (int i = 0; i < CODE_LEN; i++)
	{
		app[i] = new double[NUM_LEVEL];
	}
	double *preTx = nullptr, **rx = nullptr, **postRx = nullptr, *txFilter = nullptr, *rxFilter = nullptr;
	if (SYNCHRONOUS)
	{
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			tx[nuser] = new double[CODE_LEN];
		}
		rx = new double *[2 * CODE_LEN];
		for (int i = 0; i < 2 * CODE_LEN; i++)
		{
			rx[i] = new double[2]; // real and imaginary
		}
	}
	else
	{
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			tx[nuser] = new double[CODE_LEN*UP_RATE];
		}
		preTx = new double[CODE_LEN*UP_RATE];
		rx = new double *[CODE_LEN*UP_RATE];
		postRx = new double*[CODE_LEN*UP_RATE];
		for (int i = 0; i < CODE_LEN*UP_RATE; i++)
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
		distList = new double[2 * CODE_LEN];
		group = new int*[GROUP_SIZE];
		centroid = new double *[GROUP_SIZE];
		for (int i = 0; i < GROUP_SIZE; i++)
		{
			group[i] = new int[2 * CODE_LEN];
			centroid[i] = new double[2]; // real and imaginary
		}
		estimate = new double *[NUM_USER];
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			estimate[nuser] = new double[2]; // real and imaginary
		}
		if (EM_GMM)
		{
			softAssign = new double*[2 * CODE_LEN];
			for (int i = 0; i < 2 * CODE_LEN; i++)
			{
				softAssign[i] = new double[GROUP_SIZE];
			}
		}
	}
	double ber[SNR_NUM], bler[SNR_NUM], mse[SNR_NUM] = { 0 };
	PolarCode polar(BCT_LAYER, DATA_LEN, CRC_LEN, NUM_USER);
	FILE *result_txt = fopen("result.txt", "w");
	//---------- declaration ----------
	printf("(%d,%d) Polar-coded GDMA-BPSK system\n", CODE_LEN, DATA_LEN - NBC);
	printf("Number of users: %d\n\n", NUM_USER);
	fprintf(result_txt, "(%d,%d) Polar-coded GDMA-BPSK system\n", CODE_LEN, DATA_LEN - NBC);
	fprintf(result_txt, "Number of users: %d\n\n", NUM_USER);
	if (PCC_METHOD == 0)
	{
		printf("PCC: Bhattacharyya bound\n");
		fprintf(result_txt, "PCC: Bhattacharyya bound\n");
	}
	else if (PCC_METHOD == 1)
	{
		printf("PCC: Capacity bound\n");
		fprintf(result_txt, "PCC: Capacity bound\n");
	}
	else if (PCC_METHOD == 2)
	{
		printf("PCC: Piecewise linear approximation\n");
		fprintf(result_txt, "PCC: Piecewise linear approximation\n");
	}
	else if (PCC_METHOD == 3)
	{
		printf("PCC: Gaussain approximation\n");
		fprintf(result_txt, "PCC: Gaussain approximation\n");
	}
	if (PCC_METHOD == 4)
	{
		printf("PCC: 3GPP\n");
		fprintf(result_txt, "PCC: 3GPP\n");
	}
	else
	{
		printf("Designed SNR: %.1f dB\n", (double)DESIGHED_SNR);
		fprintf(result_txt, "Designed SNR: %.1f dB\n", (double)DESIGHED_SNR);
	}
	if (DECODING_TYPE == 1)
	{
		if (JOINT == 0)
		{
			if (LIST_SIZE == 1)
			{
				printf("Decoder: SC\n\n");
				fprintf(result_txt, "Decoder: SC\n\n");
			}
			else
			{
				printf("Decoder: SCL-%d\n\n", LIST_SIZE);
				fprintf(result_txt, "Decoder: SCL-%d\n\n", LIST_SIZE);
			}
		}
		else
		{
			if (LIST_SIZE == 1)
			{
				printf("Decoder: JPSC\n\n");
				fprintf(result_txt, "Decoder: JPSC\n\n");
			}
			else
			{
				printf("Decoder: JPSCL-%d\n\n", LIST_SIZE);
				fprintf(result_txt, "Decoder: JPSCL-%d\n\n", LIST_SIZE);
			}
		}
	}
	else
	{
		if (JOINT==0)
		{
			printf("Decoder: SPA\n\n");
			fprintf(result_txt, "Decoder: SPA\n\n");
		}
		else
		{
			printf("Decoder: J-SPA\n\n");
			fprintf(result_txt, "Decoder: J-SPA\n\n");
		}
	}
	if (CE_SCHEME == 1)
	{
		printf("Perfect CSI\n\n");
		fprintf(result_txt, "Perfect CSI\n\n");
	}
	else
	{
		if (INI_METHOD == 0)
		{
			printf("Cluster-based channel estimation:\nInitial centroids selection : LBG algo.\n");
			fprintf(result_txt, "Cluster-based channel estimation:\nInitial centroids selection : LBG algo.\n");
		}
		if (INI_METHOD == 1)
		{
			printf("Cluster-based channel estimation:\nInitial centroids selection : k-means++\n");
			fprintf(result_txt, "Cluster-based channel estimation:\nInitial centroids selection : k-means++\n");
		}
		if (INI_METHOD == 2)
		{
			printf("Cluster-based channel estimation:\nInitial centroids selection : modified k-means++\n");
			fprintf(result_txt, "Cluster-based channel estimation:\nInitial centroids selection : modified k-means++\n");
		}

		if (PROPOSAL == 1)
		{
			printf("Proposed method : 1\n");
			fprintf(result_txt, "Proposed method : 1\n");
		}
		if (PROPOSAL == 2)
		{
			printf("Proposed method : 2\n");
			fprintf(result_txt, "Proposed method : 2\n");
		}

		
		if (EM_GMM)
		{
			printf("w/ Gaussian mixture model\n");
			fprintf(result_txt, "w/ Gaussian mixture model\n");
		}
		printf("\n"); fprintf(result_txt, "\n");
	}
	if (NBC)
	{
		printf("Noncoherent block coding\n\n");
		fprintf(result_txt, "Noncoherent block coding\n\n");
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
	if (FADING_TYPE == 0)
	{
		printf("Block Fading\n");
		fprintf(result_txt, "Block Fading\n");
	}
	else if(FADING_TYPE == 1)
	{
		printf("Doppler Fading\n");
		fprintf(result_txt, "Doppler Fading\n");
	}
	else
	{
		printf("Qusic Fading\n");
		fprintf(result_txt, "Qusic Fading\n");
	}
	printf("Fading size:%d\n", FADING_SIZE);
	fprintf(result_txt, "Fading size:%d\n", FADING_SIZE);
	if (!DECODING_TYPE)
	{
		printf("Iteration size:%d\n", Iteration);
		fprintf(result_txt, "Iteration size:%d\n", Iteration);
	}
	if (!Racian)
	{
		printf("Rayleight Channel\n");
		fprintf(result_txt, "Rayleight Channel\n");
	}
	else
	{
		printf("Racian Channel ,phase = %d\n",phase);
		fprintf(result_txt, "Racian Channel,phase = %d\n",phase);	
	}
	printf("CRC_LEN = %d\n", CRC_LEN);
	fprintf(result_txt, "CRC_LEN = %d\n", CRC_LEN);
	printf("-\n\n"); fprintf(result_txt, "-\n\n");
	//---------- simulation process ----------
	polar.initialize_frozen_bits(PCC_METHOD, DESIGHED_SNR);
	for (int i = 0; i < SNR_NUM; i++)
	{
		double snrdB = SNR_START + (double)i * SNR_STEP;
		double snr = pow(10., snrdB / 10.);
		double stdDev = sqrt((double)CODE_LEN / (double)(DATA_LEN - NBC) / snr / 2.);
		long double errBlock = 0, errBit = 0;
		printf("SNR[dB] = % .1f,\n", snrdB);
		double* llr_spa = new double[2];
		llr_spa[0] = 0; llr_spa[1] = 0;
		for (int block = 1; block <= BLOCK_NUM; block++)
		{
			Transmitter(polar, data, codeword, preTx, tx, txFilter);
			EnergyProfile(chCoef);
			MultipleAccessChannel(stdDev, chCoef, tx, rx);
			if (!SYNCHRONOUS) ReceivingFilter(rx, postRx, rxFilter); 
			if (CE_SCHEME == 0) // cluster-based channel estimation
			{
				Clustering(rx, centroid, group, groupSize, distList, variation, softAssign, pow(stdDev, 2), estimate, chCoef);
				MSEComparison(chCoef, estimate, mse[i]); // user specification
			}
			MLDT(pow(stdDev, 2), chCoef, rx, app, appLlr, estimate);
			Detector(polar, data, appLlr, errBlock, errBit,app, chCoef, llr_spa);
			ber[i] = errBit / ((long double)block*(DATA_LEN - NBC)*NUM_USER);
			bler[i] = errBlock / ((long double)block*NUM_USER);
			//bler[i] = errBlock / ((long double)block * NUM_USER);
			//ber[i] = errBit / ((long double)block * NUM_USER);
			
			printf("Block# = %d, BER = %e, BLER = %e, MSE = %e\r", block, ber[i], bler[i],mse[i] / (double)(NUM_USER * block));
			//printf("Block# = %d, BLER_1 = %e, BLER_2 = %e, LLR_1 = %e, LLR_2 = %e\r", block, ber[i], bler[i], llr_spa[0]/ ((double)block * NUM_USER * DATA_LEN), llr_spa[1] / ((double)block * NUM_USER * DATA_LEN));

			if (errBlock > 1000 && block > 10000)
				break;
		}
		mse[i] /= (double)(NUM_USER*BLOCK_NUM);

		if (CE_SCHEME == 1) fprintf(result_txt, "SNR[dB] = % .1f, BER = %e, BLER = %e\n", snrdB, ber[i], bler[i]); 
		else fprintf(result_txt, "SNR[dB] = % .1f, BER = %e, BLER = %e, MSE = %e\n", snrdB, ber[i], bler[i], mse[i]); 
		printf("\n\n");
	}
	fprintf(result_txt, "\n\n");
	for (int i = 0; i < SNR_NUM; i++)
	{
		fprintf(result_txt, "%e, ", ber[i]);
	}
	fprintf(result_txt, "\n");
	for (int i = 0; i < SNR_NUM; i++)
	{
		fprintf(result_txt, "%e, ", bler[i]);
	}
	fprintf(result_txt, "\n");
	if (CE_SCHEME == 0)
	{
		for (int i = 0; i < SNR_NUM; i++)
		{
			fprintf(result_txt, "%e, ", mse[i]);
		}
		fprintf(result_txt, "\n");
	}
	fclose(result_txt);
	system("pause");
	return 0;
}