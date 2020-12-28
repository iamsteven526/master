#define _USE_MATH_DEFINES
#include <iostream>
#include "parameters.h"
#include <cstring>
#include <math.h>
using namespace std;

#pragma warning( disable : 4996 )
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
	int	**data = new int *[NUM_USER];
	int **codeword = new int *[NUM_USER];
	double **appLlr = new double *[NUM_USER];
	for (int i = 0; i < NUM_USER; i++)
	{
		data[i] = new int[DATA_LEN];
		codeword[i] = new int[CODE_LEN];
	}
	for (int i = 0; i < NUM_USER; i++)
	{
		appLlr[i] = new double[CODE_LEN];
	}
	double **app = new double *[CODE_LEN + DIFF_ENC*FFT_POINT];
	for (int i = 0; i < CODE_LEN + DIFF_ENC*FFT_POINT; i++)
	{
		app[i] = new double[NUM_LEVEL];
	}
	double ****h = new double***[NUM_USER], ****H = new double***[NUM_USER];
	for (int i = 0; i < NUM_USER; i++)
	{
		h[i] = new double**[FFT_SEGMENT + DIFF_ENC];
		H[i] = new double**[FFT_SEGMENT + DIFF_ENC];
		for (int j = 0; j < FFT_SEGMENT + DIFF_ENC; j++)
		{
			h[i][j] = new double*[FFT_POINT];
			H[i][j] = new double*[2]; // real and imaginary
			for (int k = 0; k < FFT_POINT; k++)
			{
				h[i][j][k] = new double[2]; // real and imaginary
			}
			for (int k = 0; k < 2; k++)
			{
				H[i][j][k] = new double[FFT_POINT];
			}
		}
	}
	double ***chip = new double **[NUM_USER];
	for (int i = 0; i < NUM_USER; i++)
	{
		chip[i] = new double *[FFT_SEGMENT + DIFF_ENC];
		for (int j = 0; j < FFT_SEGMENT + DIFF_ENC; j++)
		{
			chip[i][j] = new double[FFT_POINT];
		}
	}
	double ***postRx = new double **[2 * (FFT_SEGMENT + DIFF_ENC)];
	for (int i = 0; i < 2 * (FFT_SEGMENT + DIFF_ENC); i++)
	{
		postRx[i] = new double *[2]; // real and imaginary
		for (int j = 0; j < 2; j++)
		{
			if (!OVER_SAMPLE)
				postRx[i][j] = new double[FFT_POINT];
			else
				postRx[i][j] = new double[FFT_POINT * OVER_SAMPLE_RATE * UP_RATE];
		}
	}
	double ****tx = nullptr, ****preTx = nullptr, ***rx = nullptr, **tempRx = nullptr, *txFilter = nullptr, *rxFilter = nullptr;
	if (SYNCHRONOUS)
	{
		tx = new double ***[NUM_USER];
		for (int i = 0; i < NUM_USER; i++)
		{
			tx[i] = new double **[FFT_SEGMENT + DIFF_ENC];
			for (int j = 0; j < FFT_SEGMENT + DIFF_ENC; j++)
			{
				if (CP_TYPE)
				{
					tx[i][j] = new double* [CP_LEN + FFT_POINT];
					for (int k = 0; k < CP_LEN + FFT_POINT; k++)
					{
						tx[i][j][k] = new double[2]; // real and imaginary
					}
				}
				else
				{
					tx[i][j] = new double* [CP_LEN + CS_LEN + FFT_POINT];
					for (int k = 0; k < CP_LEN + CS_LEN + FFT_POINT; k++)
					{
						tx[i][j][k] = new double[2]; // real and imaginary
					}
				}
			}
		}
		rx = new double **[FFT_SEGMENT + DIFF_ENC];
		for (int i = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
		{
			if (CP_TYPE)
			{
				rx[i] = new double* [CP_LEN + FFT_POINT];
				for (int j = 0; j < CP_LEN + FFT_POINT; j++)
				{
					rx[i][j] = new double[2]; // real and imaginary
				}
			}
			else
			{
				rx[i] = new double* [CP_LEN + CS_LEN + FFT_POINT];
				for (int j = 0; j < CP_LEN + CS_LEN + FFT_POINT; j++)
				{
					rx[i][j] = new double[2]; // real and imaginary
				}
			}
		}
	}
	else
	{
		tx = new double ***[NUM_USER];
		for (int i = 0; i < NUM_USER; i++)
		{
			tx[i] = new double **[FFT_SEGMENT + DIFF_ENC];
			for (int j = 0; j < FFT_SEGMENT + DIFF_ENC; j++)
			{
				if (!OVER_SAMPLE)
				{
					if (CP_TYPE)
					{
						tx[i][j] = new double* [(CP_LEN + FFT_POINT) * UP_RATE];
						for (int k = 0; k < (CP_LEN + FFT_POINT) * UP_RATE; k++)
						{
							tx[i][j][k] = new double[2]; // real and imaginary
						}
					}
					else
					{
						tx[i][j] = new double* [(CP_LEN + CS_LEN + FFT_POINT) * UP_RATE];
						for (int k = 0; k < (CP_LEN + CS_LEN + FFT_POINT) * UP_RATE; k++)
						{
							tx[i][j][k] = new double[2]; // real and imaginary
						}
					}
				}
				else
				{
					if (CP_TYPE)
					{
						tx[i][j] = new double* [(CP_LEN + FFT_POINT) * OVER_SAMPLE_RATE * UP_RATE];
						for (int k = 0; k < (CP_LEN + FFT_POINT) * OVER_SAMPLE_RATE * UP_RATE; k++)
						{
							tx[i][j][k] = new double[2]; // real and imaginary
						}
					}
					else
					{
						tx[i][j] = new double* [(CP_LEN + CS_LEN + FFT_POINT) * OVER_SAMPLE_RATE * UP_RATE];
						for (int k = 0; k < (CP_LEN + CS_LEN + FFT_POINT) * OVER_SAMPLE_RATE * UP_RATE; k++)
						{
							tx[i][j][k] = new double[2]; // real and imaginary
						}
					}
				}
			}
		}
		rx = new double **[FFT_SEGMENT + DIFF_ENC];
		for (int i = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
		{
			if (!OVER_SAMPLE)
			{
				if (CP_TYPE)
				{
					rx[i] = new double* [(CP_LEN + FFT_POINT) * UP_RATE];
					for (int j = 0; j < (CP_LEN + FFT_POINT) * UP_RATE; j++)
					{
						rx[i][j] = new double[2]; // real and imaginary
					}
				}
				else
				{
					rx[i] = new double* [(CP_LEN + CS_LEN + FFT_POINT) * UP_RATE];
					for (int j = 0; j < (CP_LEN + CS_LEN + FFT_POINT) * UP_RATE; j++)
					{
						rx[i][j] = new double[2]; // real and imaginary
					}
				}
			}
			else
			{
				if (CP_TYPE)
				{
					rx[i] = new double* [(CP_LEN + FFT_POINT) * OVER_SAMPLE_RATE * UP_RATE];
					for (int j = 0; j < (CP_LEN + FFT_POINT) * OVER_SAMPLE_RATE * UP_RATE; j++)
					{
						rx[i][j] = new double[2]; // real and imaginary
					}
				}
				else
				{
					rx[i] = new double* [(CP_LEN + CS_LEN + FFT_POINT) * OVER_SAMPLE_RATE * UP_RATE];
					for (int j = 0; j < (CP_LEN + CS_LEN + FFT_POINT) * OVER_SAMPLE_RATE * UP_RATE; j++)
					{
						rx[i][j] = new double[2]; // real and imaginary
					}
				}
			}
		}
		preTx = new double ***[NUM_USER];
		for (int i = 0; i < NUM_USER; i++)
		{
			preTx[i] = new double **[FFT_SEGMENT + DIFF_ENC];
			for (int j = 0; j < FFT_SEGMENT + DIFF_ENC; j++)
			{
				if (!OVER_SAMPLE)
				{
					if (CP_TYPE)
					{
						preTx[i][j] = new double* [(CP_LEN + FFT_POINT) * UP_RATE];
						for (int k = 0; k < (CP_LEN + FFT_POINT) * UP_RATE; k++)
						{
							preTx[i][j][k] = new double[2]; // real and imaginary
						}
					}
					else
					{
						preTx[i][j] = new double* [(CP_LEN + CS_LEN + FFT_POINT) * UP_RATE];
						for (int k = 0; k < (CP_LEN + CS_LEN + FFT_POINT) * UP_RATE; k++)
						{
							preTx[i][j][k] = new double[2]; // real and imaginary
						}
					}
				}
				else
				{
					if (CP_TYPE)
					{
						preTx[i][j] = new double* [(CP_LEN + FFT_POINT) * OVER_SAMPLE_RATE * UP_RATE];
						for (int k = 0; k < (CP_LEN + FFT_POINT) * OVER_SAMPLE_RATE* UP_RATE; k++)
						{
							preTx[i][j][k] = new double[2]; // real and imaginary
						}
					}
					else
					{
						preTx[i][j] = new double* [(CP_LEN + CS_LEN + FFT_POINT) * OVER_SAMPLE_RATE * UP_RATE];
						for (int k = 0; k < (CP_LEN + CS_LEN + FFT_POINT) * OVER_SAMPLE_RATE * UP_RATE; k++)
						{
							preTx[i][j][k] = new double[2]; // real and imaginary
						}
					}
				}
				
			}
		}
		if (CP_TYPE)
		{
			tempRx = new double* [(CP_LEN + FFT_POINT) * UP_RATE * OVER_SAMPLE_RATE];
			for (int i = 0; i < (CP_LEN + FFT_POINT) * UP_RATE * OVER_SAMPLE_RATE; i++)
			{
				tempRx[i] = new double[2];
			}
		}
		else
		{
			tempRx = new double* [(CP_LEN + CS_LEN + FFT_POINT) * UP_RATE * OVER_SAMPLE_RATE];
			for (int i = 0; i < (CP_LEN + CS_LEN +FFT_POINT) * UP_RATE * OVER_SAMPLE_RATE; i++)
			{
				tempRx[i] = new double[2];
			}
		}
		txFilter = new double[(2 * TRUNCATION + 1)*UP_RATE];
		rxFilter = new double[(2 * TRUNCATION + 1)*UP_RATE];
		SRRCGeneration(rxFilter, 0);
	}
	double **refLlr = nullptr;
	if (DIFF_ENC)
	{
		refLlr = new double *[NUM_USER];
		for (int i = 0; i < NUM_USER; i++)
		{
			refLlr[i] = new double[FFT_POINT];
		}
	}
	int ***trellis = nullptr;
	double *preLlr = nullptr, *tempLlr = nullptr, *rxLlr = nullptr, **alpha = nullptr, **beta = nullptr, ***gamma = nullptr;
	if (TURBO_DEC)
	{
		preLlr = new double[sizeof(double)*CODE_LEN];
		tempLlr = new double[sizeof(double)*CODE_LEN];
		rxLlr = new double[sizeof(double)*CODE_LEN];
		alpha = new double*[sizeof(double)*(FFT_SEGMENT + 2)];
		beta = new double*[sizeof(double)*(FFT_SEGMENT + 2)];
		for (int i = 0; i < FFT_SEGMENT + 2; i++)
		{
			alpha[i] = new double[2];
			beta[i] = new double[2];
		}
		gamma = new double **[FFT_SEGMENT + 1];
		for (int i = 0; i < FFT_SEGMENT + 1; i++)
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
	int	*groupSize = nullptr, **group = nullptr;
	double	*variation = nullptr, *distList = nullptr, **centroid = nullptr, ****estimate = nullptr, **softAssign = nullptr;
	if (CE_SCHEME == 0)
	{
		groupSize = new int[GROUP_SIZE];
		variation = new double[GROUP_SIZE];
		distList = new double[FFT_SEGMENT];
		group = new int *[GROUP_SIZE];
		centroid = new double *[GROUP_SIZE];
		for (int i = 0; i < GROUP_SIZE; i++)
		{
			group[i] = new int[2 * (FFT_SEGMENT + DIFF_ENC)];
			centroid[i] = new double[2]; // real and imaginary
		}
		estimate = new double ***[NUM_USER];
		for (int i = 0; i < NUM_USER; i++)
		{
			estimate[i] = new double **[FFT_SEGMENT + DIFF_ENC];
			for (int j = 0; j < FFT_SEGMENT + DIFF_ENC; j++)
			{
				estimate[i][j] = new double *[FFT_POINT];
				for (int k = 0; k < FFT_POINT; k++)
				{
					estimate[i][j][k] = new double[2]; // real and imaginary
				}
			}
		}
		if (EM_GMM)
		{
			softAssign = new double*[2 * (FFT_SEGMENT + DIFF_ENC)];
			for (int i = 0; i < 2 * (FFT_SEGMENT + DIFF_ENC); i++)
			{
				softAssign[i] = new double[GROUP_SIZE];
			}
		}
	}
	double ber[SNR_NUM], bler[SNR_NUM];
	double drift[NUM_USER] = { 0 };
	int** Interleaver = new int* [NUM_USER];
	for (int i = 0; i < NUM_USER; i++)
	{
		Interleaver[i] = new int[CODE_LEN];
	}
	int error_bits_count[DATA_LEN] = { 0 };
	LDPC ldpc(LDPC_H_COL, LDPC_H_ROW);
	PolarCode polar(BCT_LAYER, DATA_LEN, CRC_LEN, NUM_USER);
	if(!CH_CODING_TYPE)
		polar.initialize_frozen_bits(PCC_METHOD, DESIGHED_SNR);
	FILE *result_txt = fopen("result.txt", "w");
	//---------- declaration ----------
	printf("OFDM-GDMA-BPSK system\n\n");
	if (CH_CODING_TYPE)
	{
		printf("(%d,%d) Polar code\n", CODE_LEN, DATA_LEN);
		fprintf(result_txt, "OFDM-GDMA-BPSK system\n\n");
		fprintf(result_txt, "(%d,%d) LDPC code\n", CODE_LEN, DATA_LEN);
		if (JCD)
		{
			printf("Joint channel decoding: G-SPA\n");
			fprintf(result_txt, "Joint channel decoding: G-SPA\n");
		}
		else
		{
			printf("Separated channel decoding: SPA\n");
			fprintf(result_txt, "Separated channel decoding: SPA\n");
		}
		if (!(DIFF_ENC && JOINT_DEC))
		{
			printf("Iteration#: %d\n\n", LDPC_IT);
			fprintf(result_txt, "Iteration#: %d\n\n", LDPC_IT);
		}
		else
		{
			printf("\n");
			fprintf(result_txt, "\n");
		}
	}
	else
	{
		printf("(%d,%d) Polar code\n", CODE_LEN, DATA_LEN);
		fprintf(result_txt, "OFDM-GDMA-BPSK system\n\n");
		fprintf(result_txt, "(%d,%d) LDPC code\n", CODE_LEN, DATA_LEN);

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
		else if (PCC_METHOD == 4)
		{
			printf("PCC: 3GPP\n");
			fprintf(result_txt, "PCC: 3GPP\n");
		}
		else if (PCC_METHOD == 5)
		{
			printf("PCC: %d -section NBC Polariztion Structure\n",FFT_POINT);
			fprintf(result_txt, "PCC: %d -section NBC Polariztion Structure\n", FFT_POINT);
		}
		else
		{
			printf("Designed SNR: %.1f dB\n", (double)DESIGHED_SNR);
			fprintf(result_txt, "Designed SNR: %.1f dB\n", (double)DESIGHED_SNR);
		}
		if (JOINT)
		{
			if (!POLAR_DECODING_TYPE)
			{
				printf("Joint channel decoding: G-SPA\n");
				fprintf(result_txt, "Joint channel decoding: G-SPA\n");
			}
			else
			{
				printf("Joint channel decoding: G-PSCL-%d\n",LIST_SIZE);
				fprintf(result_txt, "Joint channel decoding: G-PSCL-%d\n", LIST_SIZE);
			}
		}
		else
		{
			if (!POLAR_DECODING_TYPE)
			{
				printf("Separated channel decoding: SPA\n");
				fprintf(result_txt, "Separated channel decoding: SPA\n");
			}
			else
			{
				printf("Separated channel decoding: SCL-%d\n", LIST_SIZE);
				fprintf(result_txt, "Separated channel decoding: SCL-%d\n",LIST_SIZE);
			}
		}
		printf("CRC_LEN = %d", CRC_LEN);
		fprintf(result_txt, "CRC_LEN = %d", CRC_LEN);
		if (INTERLEAVER && INTERLEAVER_type)
		{
			printf(" + Interleaver\n\n");
			fprintf(result_txt, " + Interleaver\n\n");
		}
		else
		{
			printf("\n\n");
			fprintf(result_txt, "\n\n");
		}
	}
	if (CH_TYPE == 2)
	{
		printf("Time-varing channel\n");
		printf("UE speed: %d km/hour\n\n", UE_SPEED);
		fprintf(result_txt, "Time-varing channel\n");
		fprintf(result_txt, "UE speed: %d km/hour\n\n", UE_SPEED);
	}
	printf("FFT point: %d, CP length: %d, channel taps: %d\n", FFT_POINT, CP_LEN, TAP_NUM);
	printf("Number of users: %d\n\n", NUM_USER);
	fprintf(result_txt, "FFT point: %d, CP length: %d, channel tap: %d\n", FFT_POINT, CP_LEN, TAP_NUM);
	fprintf(result_txt, "Number of users: %d\n\n", NUM_USER);
	if (CE_SCHEME == 1)
	{
		printf("Perfect CSI\n\n");
		fprintf(result_txt, "Perfect CSI\n\n");
		if (NBC)
		{
			cout << "perfect csi + NBC error !!\n\n";
		}
	}
	else
	{
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
		if (SLIDING)
		{
			printf("Sliding window with size: %d\n", WINDOW_SIZE);
			fprintf(result_txt, "Sliding window with size: %d\n", WINDOW_SIZE);
		}
		
		printf("\n"); fprintf(result_txt, "\n");
	}
	
	if (CH_CODING_TYPE && DIFF_ENC)
	{
		printf("Differential encoding,\n");
		fprintf(result_txt, "Differential encoding,\n");
	}
	if (CH_CODING_TYPE && DIFF_ENC && DIFF_DEC)
	{
		printf("w/ differential decoding\n\n");
		fprintf(result_txt, "w/ differential decoding\n\n");
	}
	if (CH_CODING_TYPE && DIFF_ENC && TURBO_DEC)
	{
		printf("w/ turbo processing\n");
		printf("Turbo iterations: %d\n\n", TURBO_IT);
		fprintf(result_txt, "w/ turbo processing\n");
		fprintf(result_txt, "Turbo iterations: %d\n\n", TURBO_IT);
	}
	if (CH_CODING_TYPE && DIFF_ENC && JOINT_DEC)
	{
		printf("w/ joint decoding\n");
		printf("Joint/inner iterations: %d/%d\n\n", JOINT_IT, INNER_IT);
		fprintf(result_txt, "w/ joint decoding\n");
		fprintf(result_txt, "Joint/inner iterations: %d/%d\n\n", JOINT_IT, INNER_IT);
	}
	if (SYNCHRONOUS)
	{
		printf("Synchronous system\n");
		fprintf(result_txt, "Synchronous system\n");
	}
	else
	{
		printf("Asynchronous system\n");
		printf("Drifting range: [%.3f, %.3f]\n", -0.5*DRIFT_RANGE, 0.5*DRIFT_RANGE);
		fprintf(result_txt, "Asynchronous system\n");
		fprintf(result_txt, "Drifting range: [%.3f, %.3f]\n", -0.5*DRIFT_RANGE, 0.5*DRIFT_RANGE);
	}
	printf("Channel, K=%d\n",K);
	fprintf(result_txt, "Channel, K=%d\n", K);

	
	if (CP_TYPE)
	{
		printf("OFDM : CP\n");
		fprintf(result_txt, "OFDM : CP\n");
	}
	else
	{
		printf("OFDM : CP + CS\n");
		fprintf(result_txt, "OFDM : CP + CS\n");
	}

	if(OVER_SAMPLE)
	{
		if(PULSE_SHAPE)
		{
			printf("Tx : Over Sampling:%d ; Pulse shaping:%d\n", OVER_SAMPLE_RATE,UP_RATE);
			fprintf(result_txt, "Tx : Over Sampling:%d ; Pulse shaping:%d\n", OVER_SAMPLE_RATE, UP_RATE);
		}
		else
		{
			printf("Tx : Over Sampling\n");
			fprintf(result_txt, "Tx : Over Sampling\n");
		}
	}
	else
	{
		if (PULSE_SHAPE)
		{
			printf("Tx : Pulse shaping:%d\n", UP_RATE);
			fprintf(result_txt, "Tx : Pulse shaping:%d\n", UP_RATE);
		}
		else
		{
			if (!SYNCHRONOUS)
			{
				cout << "Setting Error !!";
				system("pause");
			}
		}
	}
	printf("-\n"); fprintf(result_txt, "-\n");
	//---------- simulation process ---------
	for (int i = 0; i < SNR_NUM; i++)
	{
		double snrdB = SNR_START + (double)i * SNR_STEP;
		double snr = pow(10., snrdB / 10.);
		double stdDev = sqrt((0.5 / snr) * ((double)(CODE_LEN + DIFF_ENC*FFT_POINT) / (DATA_LEN - (NBC * (FFT_POINT-1)))) * ((double)(CP_LEN + FFT_POINT) / FFT_POINT));
		long double errCount[2] = { 0 }; // 0: bit error; 1: block error
		printf("SNR[dB] = % .1f,\n", snrdB);
		
		double llr = 0;
		for (int block = 1; block <= BLOCK_NUM; block++)
		{
			EnergyProfile(h, H);
			Encoder(ldpc, polar, data, codeword, Interleaver);
			Modulator(codeword, chip);
			MultiCarrierMapper(chip, tx);
			if (!SYNCHRONOUS)
			{
				UpSampling(tx, preTx);
				if (PULSE_SHAPE)
					PulseShaping(preTx, tx, txFilter, drift);
			}
			MultipleAccessChannel(stdDev, h, tx, rx, drift, H);
			if (!SYNCHRONOUS && PULSE_SHAPE) ReceivingFilter(rx, tempRx, rxFilter);
			MultiCarrierDemapper(rx, postRx, drift);

			if (CE_SCHEME == 0)
			{
				Clustering(postRx, centroid, group, groupSize, distList, variation, softAssign, pow(stdDev, 2), estimate);
				UserIdentification(H, estimate);
			}
			MLDT(ldpc, pow(stdDev, 2), H, postRx, app, appLlr, refLlr, estimate);
			//if (CH_CODING_TYPE && DIFF_ENC && TURBO_DEC) TurboProcessor(ldpc, trellis, tempLlr, refLlr, rxLlr, preLlr, appLlr, alpha, beta, gamma);
			Detector(ldpc, polar, data, appLlr, refLlr, errCount, app, Interleaver, error_bits_count);
			ber[i] = NBC ? errCount[0] / ((long double)block * (DATA_LEN - (FFT_POINT - 1)) * NUM_USER) : errCount[0] / ((long double)block * DATA_LEN * NUM_USER);
			bler[i] = errCount[1] / ((long double)block * NUM_USER);
			printf("Block# = %d, BER = %e, BLER = %e\r", block, ber[i], bler[i]);

			if (errCount[1] > 1000 && block > 2000)
			{
				/*for (int i = 0; i < DATA_LEN; i++)
					cout << error_bits_count[i] << endl;
				system("pause");*/
				break;
			}
		}
		fprintf(result_txt, "SNR[dB] = % .1f, BER = %e, BLER = %e\n", snrdB, ber[i], bler[i]);
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
	fclose(result_txt);
	system("pause");
	return 0;
}