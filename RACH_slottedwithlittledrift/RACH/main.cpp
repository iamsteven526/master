#define _USE_MATH_DEFINES
#include <iostream>
#include "parameters.h"
#include <cstring>
#include <cmath>
using namespace std;

#pragma warning( disable : 4996 )
int main()
{
	//---------- check for overflow ----------
	if ((long double)BLOCK_NUM*DATA_LEN*NUM_USER < 0)
	{
		cout << "OVERFLOW" << endl;
		//system("pause");
		return 0;
	}
	//---------- memory allocation
	int	***data = new int **[TEST_FRAME_SIZE];
	int ***codeword = new int **[TEST_FRAME_SIZE];
	double ***appLlr = new double **[TEST_FRAME_SIZE];
	for (int t = 0; t < TEST_FRAME_SIZE; t++)
	{
		data[t] = new int* [NUM_USER];
		codeword[t] = new int* [NUM_USER];
		appLlr[t] = new double* [NUM_USER];
		for (int i = 0; i < NUM_USER; i++)
		{
			data[t][i] = new int[DATA_LEN];
			codeword[t][i] = new int[CODE_LEN];
		}
		for (int i = 0; i < NUM_USER; i++)
		{
			appLlr[t][i] = new double[CODE_LEN];
		}
	}
	double**** app = new double*** [TEST_FRAME_SIZE];
	for (int t = 0; t < TEST_FRAME_SIZE; t++)
	{
		app[t] = new double** [NUM_USER];
		for(int nuser=0; nuser<NUM_USER; nuser++)
		{
			app[t][nuser] = new double* [CODE_LEN + DIFF_ENC * FFT_POINT];
			for (int i = 0; i < CODE_LEN + DIFF_ENC * FFT_POINT; i++)
			{
				app[t][nuser][i] = new double[NUM_LEVEL];
			}
		}
	}
	double***** h = new double**** [TEST_FRAME_SIZE], ***** H = new double**** [TEST_FRAME_SIZE];
	for (int t = 0; t < TEST_FRAME_SIZE; t++)
	{
		h[t] = new double*** [NUM_USER];
		H[t] = new double*** [NUM_USER];
		for (int i = 0; i < NUM_USER; i++)
		{
			h[t][i] = new double** [FFT_SEGMENT + DIFF_ENC];
			H[t][i] = new double** [FFT_SEGMENT + DIFF_ENC];
			for (int j = 0; j < FFT_SEGMENT + DIFF_ENC; j++)
			{
				h[t][i][j] = new double* [FFT_POINT];
				H[t][i][j] = new double* [2]; // real and imaginary
				for (int k = 0; k < FFT_POINT; k++)
				{
					h[t][i][j][k] = new double[2]; // real and imaginary
				}
				for (int k = 0; k < 2; k++)
				{
					H[t][i][j][k] = new double[FFT_POINT];
				}
			}
		}
	}

	double*** h_preamble = new double** [NUM_USER];
	double*** H_preamble = new double** [NUM_USER];
	for (int i = 0; i < NUM_USER; i++)
	{
		h_preamble[i] = new double* [PREABLE_LEN];
		H_preamble[i] = new double* [2];
		for (int j = 0; j < PREABLE_LEN; j++)
		{
			h_preamble[i][j] = new double[2];
		}
		for (int j = 0; j < 2; j++)
		{
			H_preamble[i][j] = new double[PREABLE_LEN];
		}
	}
	
	double**** chip = new double*** [TEST_FRAME_SIZE];
	for (int t = 0; t < TEST_FRAME_SIZE; t++)
	{
		chip[t] = new double** [NUM_USER];
		for (int i = 0; i < NUM_USER; i++)
		{
			chip[t][i] = new double* [FFT_SEGMENT + DIFF_ENC];
			for (int j = 0; j < FFT_SEGMENT + DIFF_ENC; j++)
			{
				chip[t][i][j] = new double[FFT_POINT];
			}
		}
	}
	double***** postRx = new double**** [TEST_FRAME_SIZE];
	for (int t = 0; t < TEST_FRAME_SIZE; t++)
	{
		postRx[t] = new double*** [NUM_USER];
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			postRx[t][nuser] = new double** [FFT_SEGMENT + DIFF_ENC];
			for (int i = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
			{
				postRx[t][nuser][i] = new double* [2]; // real and imaginary
				for (int j = 0; j < 2; j++)
				{
					if (!OVER_SAMPLE)
						postRx[t][nuser][i][j] = new double[FFT_POINT];
					else
						postRx[t][nuser][i][j] = new double[FFT_POINT * OVER_SAMPLE_RATE * UP_RATE];
				}
			}
		}
	}
	
	double *****tx = nullptr, *****preTx = nullptr, **rx = nullptr, ***tempRx = nullptr, *txFilter = nullptr, *rxFilter = nullptr, *****ClassifyRx = nullptr;
	
	tx = new double**** [TEST_FRAME_SIZE];
	ClassifyRx = new double**** [TEST_FRAME_SIZE];
	for (int t = 0; t < TEST_FRAME_SIZE; t++)
	{
		tx[t] = new double*** [NUM_USER];
		ClassifyRx[t] = new double*** [NUM_USER];

		for (int i = 0; i < NUM_USER; i++)
		{
			tx[t][i] = new double** [FFT_SEGMENT + DIFF_ENC];
			ClassifyRx[t][i] = new double** [FFT_SEGMENT + DIFF_ENC];
			for (int j = 0; j < FFT_SEGMENT + DIFF_ENC; j++)
			{
				if (CP_TYPE)
				{
					tx[t][i][j] = new double* [CP_LEN + FFT_POINT];
					ClassifyRx[t][i][j] = new double* [CP_LEN + FFT_POINT];
					for (int k = 0; k < CP_LEN + FFT_POINT; k++)
					{
						tx[t][i][j][k] = new double[2]; // real and imaginary
						ClassifyRx[t][i][j][k] = new double[2]; // real and imaginary
					}
				}
				else
				{
					tx[t][i][j] = new double* [CP_LEN + CS_LEN + FFT_POINT];
					ClassifyRx[t][i][j] = new double* [CP_LEN + CS_LEN + FFT_POINT];
					for (int k = 0; k < CP_LEN + CS_LEN + FFT_POINT; k++)
					{
						tx[t][i][j][k] = new double[2]; // real and imaginary
						ClassifyRx[t][i][j][k] = new double[2]; // real and imaginary
					}
				}
			}
		}
	}

	if (CP_TYPE)
	{
		rx = new double* [((FFT_SEGMENT + DIFF_ENC) * (CP_LEN + FFT_POINT) + (CP_LEN + CS_LEN + PREABLE_LEN) * PREABLE) * TEST_FRAME_SIZE * Unit];
		for (int j = 0; j < ((FFT_SEGMENT + DIFF_ENC) * (CP_LEN + FFT_POINT) + (CP_LEN + CS_LEN + PREABLE_LEN) * PREABLE)* TEST_FRAME_SIZE * Unit; j++)
		{
			rx[j] = new double[2]; // real and imaginary
		}
	}
	else
	{
		rx = new double* [((FFT_SEGMENT + DIFF_ENC) * (CP_LEN + CS_LEN + FFT_POINT) + (CP_LEN + CS_LEN + PREABLE_LEN) * PREABLE) * TEST_FRAME_SIZE * Unit];
		for (int j = 0; j < ((FFT_SEGMENT + DIFF_ENC) * (CP_LEN + CS_LEN + FFT_POINT) + (CP_LEN + CS_LEN + PREABLE_LEN) * PREABLE) * TEST_FRAME_SIZE * Unit; j++)
		{
			rx[j] = new double[2]; // real and imaginary
		}
	}
	if(PULSE_SHAPE)
	{
		txFilter = new double[(2 * TRUNCATION + 1) * UP_RATE];
		rxFilter = new double[(2 * TRUNCATION + 1) * UP_RATE];
		SRRCGeneration(rxFilter, 0);
		SRRCGeneration(txFilter, 0);
	}

	double ***refLlr = nullptr;
	if (DIFF_ENC)
	{
		refLlr = new double** [TEST_FRAME_SIZE];
		for (int t = 0; t < TEST_FRAME_SIZE; t++)
		{
			refLlr[t] = new double* [NUM_USER];
			for (int i = 0; i < NUM_USER; i++)
			{
				refLlr[t][i] = new double[FFT_POINT];
			}
		}
	}
	
	int	*groupSize = nullptr, **group = nullptr;
	double	*variation = nullptr, *distList = nullptr, **centroid = nullptr, *****estimate = nullptr, **cluster_sample = nullptr, ** softAssign = nullptr;
	if (CE_SCHEME == 0)
	{
		groupSize = new int[GROUP_SIZE];
		variation = new double[GROUP_SIZE];
		distList = new double[2 * FFT_SEGMENT];
		group = new int *[GROUP_SIZE];
		centroid = new double *[GROUP_SIZE];
		for (int i = 0; i < GROUP_SIZE; i++)
		{
			group[i] = new int[FFT_SEGMENT + DIFF_ENC];
			centroid[i] = new double[2]; // real and imaginary
		}
		estimate = new double ****[TEST_FRAME_SIZE];
		for (int t = 0; t < TEST_FRAME_SIZE; t++)
		{
			estimate[t] = new double ***[NUM_USER];
			for (int i = 0; i < NUM_USER; i++)
			{
				estimate[t][i] = new double **[FFT_SEGMENT + DIFF_ENC];
				for (int j = 0; j < FFT_SEGMENT + DIFF_ENC; j++)
				{
					estimate[t][i][j] = new double *[FFT_POINT];
					for (int k = 0; k < FFT_POINT; k++)
					{
						estimate[t][i][j][k] = new double[2]; // real and imaginary
					}
				}
			}
		}

		cluster_sample = new double* [2 * (FFT_SEGMENT + DIFF_ENC)];
		for (int i = 0; i < 2 * (FFT_SEGMENT + DIFF_ENC); i++)
		{
			cluster_sample[i] = new double[2];
		}

		if (EM_GMM)
		{
			softAssign = new double* [2 * FFT_SEGMENT];
			for (int i = 0; i < 2 * FFT_SEGMENT; i++)
			{
				softAssign[i] = new double[GROUP_SIZE];
			}
		}
	}
	double S[G_NUM], mse[G_NUM], Time_deviation[G_NUM], Time_error[G_NUM];
	double drift[NUM_USER] = { 0 };
	int *packet_num = new int[NUM_USER], *estimate_packet_num = new int[NUM_USER];
	int** packet_time = new int* [NUM_USER];
	int* Time_Offset_Error = new int [(CP_LEN + CS_LEN) * Unit];
	bool **estimate_packet_time = new bool* [NUM_USER];
	for (int i = 0; i < NUM_USER; i++)
	{
		packet_time[i] = new int[TEST_FRAME_SIZE];
		estimate_packet_time[i] = new bool[TEST_FRAME_SIZE];
	}
	int*** Cluster_num = new int** [TEST_FRAME_SIZE], **** Cluster_gain = new int*** [TEST_FRAME_SIZE];
	for (int t = 0; t < TEST_FRAME_SIZE; t++)
	{
		Cluster_num[t] = new int* [NUM_USER];
		Cluster_gain[t] = new int** [NUM_USER];
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			Cluster_num[t][nuser] = new int [FFT_SEGMENT + DIFF_ENC];
			Cluster_gain[t][nuser] = new int* [FFT_SEGMENT + DIFF_ENC];
			for (int s = 0; s < FFT_SEGMENT + DIFF_ENC; s++)
			{
				Cluster_gain[t][nuser][s] = new int [NUM_USER];
			}
		}
	}
	double*** ZC_sequence = new double** [Root_Num];
	for (int i = 0; i < Root_Num; i++)
	{
		ZC_sequence[i] = new double* [PREABLE_LEN];
		for (int j = 0; j < PREABLE_LEN; j++)
		{
			ZC_sequence[i][j] = new double[2];
		}
	}


	LDPC ldpc(LDPC_H_COL, LDPC_H_ROW);
	FILE *result_txt = fopen("result.txt", "w");
	//---------- declaration ----------
	printf("RACH-OFDM-GDMA-BPSK system\n\n");
	printf("(%d,%d) LDPC code\n", CODE_LEN, DATA_LEN);
	fprintf(result_txt, "RACH-OFDM-GDMA-BPSK system\n\n");
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
	if (CH_TYPE == 2)
	{
		printf("Time-varing channel\n");
		printf("UE speed: %d km/hour\n\n", UE_SPEED);
		fprintf(result_txt, "Time-varing channel\n");
		fprintf(result_txt, "UE speed: %d km/hour\n\n", UE_SPEED);
	}
	printf("FFT point: %d, CP length: %d, channel taps: %d\n", FFT_POINT, CP_LEN, TAP_NUM);
	printf("Max Number of users: %d\n\n", NUM_USER);
	fprintf(result_txt, "FFT point: %d, CP length: %d, channel tap: %d\n", FFT_POINT, CP_LEN, TAP_NUM);
	fprintf(result_txt, "Number of users: %d\n\n", NUM_USER);
	if (!Time_Estimate)
	{
		printf("Perfect Time\n\n");
		fprintf(result_txt, "Perfect Time\n\n");
	}
	else
	{
		printf("Estimated Time\n\n");
		fprintf(result_txt, "Estimated Time\n\n");
	}
	
	if (CE_SCHEME == 1)
	{
		printf("Perfect CSI\n\n");
		fprintf(result_txt, "Perfect CSI\n\n");
	}
	else
	{
		if (CE_METHOD == 0)
		{
			printf("Preamble Estimation Only\n");
			fprintf(result_txt, "Preamble Estimation Only\n");
		}
		if (CE_METHOD == 1)
		{
			printf("Preamble + GMM\n");
			fprintf(result_txt, "Preamble + GMM\n");
		}
		
		if (SLIDING)
		{
			printf("Sliding window with size: %d\n", WINDOW_SIZE);
			fprintf(result_txt, "Sliding window with size: %d\n", WINDOW_SIZE);
		}
		printf("\n"); fprintf(result_txt, "\n");
	}
	if (DIFF_ENC)
	{
		printf("Differential encoding,\n");
		fprintf(result_txt, "Differential encoding,\n");
	}
	if (DIFF_ENC && DIFF_DEC)
	{
		printf("w/ differential decoding\n\n");
		fprintf(result_txt, "w/ differential decoding\n\n");
	}
	if (DIFF_ENC && TURBO_DEC)
	{
		printf("w/ turbo processing\n");
		printf("Turbo iterations: %d\n\n", TURBO_IT);
		fprintf(result_txt, "w/ turbo processing\n");
		fprintf(result_txt, "Turbo iterations: %d\n\n", TURBO_IT);
	}
	if (DIFF_ENC && JOINT_DEC)
	{
		printf("w/ joint decoding\n");
		printf("Joint/inner iterations: %d/%d\n\n", JOINT_IT, INNER_IT);
		fprintf(result_txt, "w/ joint decoding\n");
		fprintf(result_txt, "Joint/inner iterations: %d/%d\n\n", JOINT_IT, INNER_IT);
	}


	if (ALOHA)
	{
		printf("ALOHA + MLDT, ");
		fprintf(result_txt, "ALOHA + MLDT, \t");
	}
	else
	{
		printf("ALOHA, ");
		fprintf(result_txt, "ALOHA, ");
	}
	if (SLOTED)
	{
		printf("Slotted system, Test Frame : %d\n", TEST_FRAME_SIZE);
		fprintf(result_txt, "SLOTED systeme\n");
	}
	else
	{
		printf("NON Slotted system, Test Frame : %d\n", TEST_FRAME_SIZE);
		fprintf(result_txt, "NON SLOTED system\n");
	}

	printf("Channel, K=%d\n",K);
	fprintf(result_txt, "Channel, K=%d\n", K);

	
	if (CP_TYPE)
	{
		printf("OFDM : CP, CP Len : %d\n", CP_LEN);
		fprintf(result_txt, "OFDM : CP, CP Len : %d\n", CP_LEN);
	}
	else
	{
		printf("OFDM : CP + CS, CP Len : %d, CS Len : %d\n", CP_LEN, CS_LEN);
		fprintf(result_txt, "OFDM : CP + CS, CP Len : %d, CS Len : %d\n", CP_LEN, CS_LEN);
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
	}

	if (PREABLE)
	{
		printf("Preamble Len: %d\n", PREABLE_LEN);
		fprintf(result_txt, "Preamble Len: %d\n", PREABLE_LEN);
	}
	printf("-\n"); fprintf(result_txt, "-\n");
	//---------- simulation process ----------
	ZC_Generater(ZC_sequence);
	for (int i = 0; i < G_NUM; i++)
	{
		double G = G_START + (double)i * G_STEP;
		double snrdB = 30;
		//double G = 0.1;
		//double snrdB = 20 + (double)i * 2;

		double snr = pow(10., snrdB / 10.);
		double stdDev = sqrt((0.5 / snr) * ((double)(CODE_LEN + DIFF_ENC*FFT_POINT) / DATA_LEN) * ((double)(CP_LEN + FFT_POINT) / FFT_POINT));
		long double errCount[2] = { 0 }; // 0: bit error; 1: block error
		long double packet_sum = 0;
		int pre_time = 0;
		
		printf("G = % .1f / packet_dur; SNR = % .1f dB\n", G, snrdB);

		double** segment_num=new double*[2], ** segment_error=new double*[2];
		/*for (int k = 0; k < 2; k++)
		{
			segment_num[k] = new double[20*62];
			segment_error[k] = new double[20 * 62];
		}
		for (int k = 0; k < 2; k++)
		{
			for (int kk = 0; kk < 20 * 62; kk++)
			{
				segment_num[k][kk] = 0;
				segment_error[k][kk] = 0;
			}
		}*/

		/*for (int i = 0; i < 40; i++)
		{
			Time_Offset_Error[i] = 0;
		}*/
		for (int block = 1; block <= BLOCK_NUM; block++)
		{
			Packet_generater(packet_num, packet_time, packet_sum,G, errCount);
			
			EnergyProfile(h, H, packet_num, packet_time, h_preamble, H_preamble);
			
			Encoder(ldpc, data, codeword, packet_num);
			Modulator(codeword, chip, packet_num);
			MultiCarrierMapper(chip, tx, packet_num);
			
			MultipleAccessChannel(stdDev, h, tx, rx, drift, H, packet_num, packet_time, txFilter, ZC_sequence, h_preamble);
		
			Time_Num_Estimation(rx, estimate_packet_time,Time_deviation[i], Time_error[i], packet_num, packet_time, Time_Offset_Error, estimate, ZC_sequence);
			
			
			ReceivingClassify(rx, ClassifyRx, packet_num, packet_time, Cluster_num, Cluster_gain);
			
			MultiCarrierDemapper(ClassifyRx, postRx, packet_num);
			
			
			if (CE_SCHEME == 0)
			{
				Clustering(postRx, pow(stdDev, 2), estimate, H, Cluster_num, packet_num,mse[i], estimate_packet_time, cluster_sample, centroid, group, groupSize, variation, packet_time, softAssign);
				//UserIdentification(H, estimate, mse[i]);
			}
			
			
			MLDT(ldpc, stdDev, H, postRx, app, appLlr, refLlr, estimate, packet_num, Cluster_num, Cluster_gain, snrdB);

			Detector(ldpc, data, appLlr, refLlr, errCount, packet_num, estimate_packet_time, segment_num, segment_error, packet_time);
			
			
			S[i] = (packet_sum - errCount[1]) / ((long double) TEST_FRAME_SIZE * block);
			//S[i] = (packet_sum - errCount[1]) / packet_sum;
			//S[i] = errCount[1] / block;
			
			printf("Block# = %d, S = %e, G = %e, mse = %e, T_D = %e, T_E = %e \r", block, S[i], G, mse[i] / block, Time_deviation[i] / (packet_sum - Time_error[i]), Time_error[i]/packet_sum);
			//cout << endl;
			//system("pause");
			//printf("Block# = %d, E0 = %d, E1 = %d, E2 = %d, E3 = %d \r", block, Time_Offset_Error[0], Time_Offset_Error[1], Time_Offset_Error[2], Time_Offset_Error[3]);
			
			/*if (block == 100000)
			{
				for (int i = 0; i < 40; i++)
				{
					cout << Time_Offset_Error[i] << " ";
				}
				system("pause");
				for (int k = 0; k < 20 * 62; k++)
				{
					if (segment_num[1][k] == 0)
					{
						cout << 0 << " ";
						//fprintf(result_txt, "%d ", 0);
					}
					else
					{
						cout <<  segment_error[1][k]/segment_num[1][k] << " ";
						//fprintf(result_txt, "%e ", segment_error[1][k] / segment_num[1][k]);
					}
				}
			}*/
			if (block > 10000)
				break;

			//if (errCount[1] > 1000 || block > 1000)
			//	break;

			
		}
		
		/*for (int k = 0; k < 20 * 62; k++)
		{
			if (segment_num[1][k] == 0)
			{
				cout << 0 << " ";
				fprintf(result_txt, "%d ", 0);
			}
			else
			{
				cout <<  segment_error[1][k]/segment_num[1][k] << " ";
				//fprintf(result_txt, "%e ", segment_error[1][k] / segment_num[1][k]);
			}
		}*/
		//fclose(result_txt);
		//system("pause");
		mse[i] /= BLOCK_NUM;
		fprintf(result_txt, "Packet Rate = % .1f, BER = %e, BLER = %e, MSE = %e, T_D = %e, T_E = %e \n", S[i], G, mse[i] , Time_deviation[i] / (packet_sum - Time_error[i]), Time_error[i] / packet_sum);
		printf("\n\n");
	}
	fprintf(result_txt, "\n\n");
	for (int i = 0; i < G_NUM; i++)
	{
		fprintf(result_txt, "%e, ", S[i]);
	}
	fprintf(result_txt, "\n");
	for (int i = 0; i < G_NUM; i++)
	{
		fprintf(result_txt, "%e, ", mse[i]);
	fprintf(result_txt, "\n");
	fclose(result_txt);
	//system("pause");
	return 0;
}
	}