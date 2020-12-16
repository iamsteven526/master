#include <iostream>
#include "parameters.h"
#include <cmath>
using namespace std;
#pragma warning (disable : 4996)

//--------- scma matrix-------------//can be generalize
int		 scma_matrix[SCMA_SOURCE][SCMA_SOURCE * NUM_USER / SCMA_USER_SOURCE] = { {0,1,1,0,1,0}, {1,0,1,0,0,1}, {0,1,0,1,0,1}, {1,0,0,1,1,0} };
int		 coef_idx[SCMA_SOURCE][SCMA_SOURCE * NUM_USER / SCMA_USER_SOURCE] = { {-1,0,0,-1,0,-1}, {0,-1,1,-1,-1,0}, {-1,1,-1,0,-1,1}, {1,-1,-1,1,1,-1} };

//12users BPSK 2020/02/19

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
	int** data = new int* [NUM_USER * SCMA_SOURCE / SCMA_USER_SOURCE];
	double** tx = new double* [NUM_USER * SCMA_SOURCE / SCMA_USER_SOURCE];
	double** appLlr = new double* [NUM_USER * SCMA_SOURCE / SCMA_USER_SOURCE];
	double*** chCoef = new double** [2];//BPSK user1
	double*** chCoef2 = new double** [2];//BPSK user2

    for (int resource = 0; resource < 2; resource++){
        chCoef[resource] = new double* [NUM_USER * SCMA_SOURCE / SCMA_USER_SOURCE];//BPSK user1
	    chCoef2[resource] = new double* [NUM_USER * SCMA_SOURCE / SCMA_USER_SOURCE];//BPSK user2	
		for (int nuser = 0; nuser < NUM_USER * SCMA_SOURCE / SCMA_USER_SOURCE; nuser++)
	    {
		    chCoef[resource][nuser] = new double[5];
		    chCoef2[resource][nuser] = new double[5];
		}	
	}


	for (int nuser = 0; nuser < NUM_USER * SCMA_SOURCE / SCMA_USER_SOURCE; nuser++)
	{
		data[nuser] = new int[BLOCK_LEN * Qm];
		tx[nuser] = new double[BLOCK_LEN * Qm];
		appLlr[nuser] = new double[BLOCK_LEN * Qm];

	}
	double** app = new double* [SCMA_SOURCE];//TODO ?? think why?? 20200809
	for (int i = 0; i < SCMA_SOURCE; i++)
	{
		app[i] = new double[NUM_LEVEL];
	}
	cout << NUM_LEVEL;
	double*** rx = new double** [SCMA_SOURCE];
	for (int j = 0; j < SCMA_SOURCE; ++j) {
		rx[j] = new double* [2 * BLOCK_LEN];
		for (int i = 0; i < (2* BLOCK_LEN); i++)
		{
			rx[j][i] = new double[2]; // real and imaginary
		}
	}

	int* groupSize = nullptr, ** group = nullptr;
	double *variation = nullptr, *distList = nullptr, **centroid = nullptr, **estimate = nullptr, **softAssign = nullptr, ***finalestimate = nullptr, ***finalestimate2 = nullptr;
	
	if (CE_METHOD == 0)
	{
		groupSize = new int[GROUP_SIZE];
		variation = new double[GROUP_SIZE];
		distList = new double[2*BLOCK_LEN];
		group = new int* [GROUP_SIZE];
		centroid = new double* [GROUP_SIZE];
		for (int i = 0; i < GROUP_SIZE; i++)
		{
			group[i] = new int[BLOCK_LEN];
			centroid[i] = new double[2]; // real and imaginary
		}
		estimate = new double* [NUM_USER* SCMA_SOURCE / SCMA_USER_SOURCE];
        finalestimate = new double** [2];
		finalestimate2 = new double** [2];
        for (int resource = 0; resource < 2; resource++){
		    finalestimate[resource] = new double* [NUM_USER* SCMA_SOURCE / SCMA_USER_SOURCE];
		    finalestimate2[resource] = new double* [NUM_USER* SCMA_SOURCE / SCMA_USER_SOURCE];
			for (int i = 0; i < (NUM_USER * SCMA_SOURCE / SCMA_USER_SOURCE); i++)
		    {
			    finalestimate[resource][i] = new double[2]; // real and imaginary
			    finalestimate2[resource][i] = new double[2]; // real and imaginary
			}
		}



		for (int i = 0; i < (NUM_USER * SCMA_SOURCE / SCMA_USER_SOURCE); i++)
		{
			estimate[i] = new double[2]; // real and imaginary
		}
		if (EM_GMM)
		{
			softAssign = new double* [2 * BLOCK_LEN];
			for (int i = 0; i < 2 * BLOCK_LEN; i++)
			{
				softAssign[i] = new double[GROUP_SIZE];
			}
		}
	}
	int* known_drift = nullptr;
	//if (COLLISION == 0)
	known_drift = new int[2 * NUM_USER * SCMA_SOURCE / SCMA_USER_SOURCE];

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
	double ber[SNR_NUM], mse[SNR_NUM] = { 0 }, mse_centroid[SNR_NUM] = { 0 };
	long double itCount[SNR_NUM] = { 0 };


	FILE* result_txt = fopen("result.txt", "w");
	//---------- declaration ----------
	printf("Uncoded GDMA-BPSK system\n");
	printf("Number of users: %d\n\n", NUM_USER);
	printf("Block length: %d\n\n", BLOCK_LEN);
	printf("Cluster length: %d\n\n", Cluster_Len);
	fprintf(result_txt, "GDMA-BPSK system\n");
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
		double variance = 0.5 / snr;
		double stdDev = sqrt(0.25 / snr);
		long double error = 0;
		printf("SNR[dB] = %.1f\n", snrdB);
		for (int block = 1; block <= BLOCK_NUM; ++block)
		{
			Transmitter(data, tx);
			EnergyProfile(chCoef); 
			EnergyProfile(chCoef2);
			for (int source_id = 0; source_id < SCMA_SOURCE; ++source_id) {//?				
				MultipleAccessChannel(stdDev, chCoef, chCoef2, tx, rx[source_id], source_id, scma_matrix, coef_idx);
			}
			//if (!SYNCHRONOUS) ReceivingFilter(rx, postRx, rxFilter, txFilter, pilot,chCoef);
			if (CE_METHOD == 0) // cluster-based channel estimation
			{
				//KmeansClustering(rx, centroid, group, groupSize, distList, variation, softAssign, variance, estimate, itCount[i], chCoef);
				for (int source_id = 0; source_id < SCMA_SOURCE; ++source_id) {
					//cout << block << ',' << source_id << endl;
				    Clustering(rx[source_id], centroid, group, groupSize, distList, variation, softAssign, variance, estimate, itCount[i],chCoef,known_drift);
				    //debug

					//for(int qq = 0;qq<6;++qq){
				     //   cout << "cluster  " << 1.414*estimate[qq][0] << " " << 1.414*estimate[qq][1] << endl;
			        //}
					
					SCMA_MSEComparison(chCoef,chCoef2, estimate,finalestimate,finalestimate2,scma_matrix[source_id],coef_idx[source_id]);
					//MSEComparison(chCoef, estimate, mse[i], rx[source_id], centroid,known_drift); // user specification
				}
				//CentroidMSEComparison(chCoef, estimate, mse_centroid[i], rx, centroid);
			}
			//debug
			//for(int qq = 0;qq<6;++qq){
			//	cout << chCoef[0][qq][0] << " " << chCoef[0][qq][1]<< "         2bb       " << chCoef2[0][qq][0] << " " << chCoef2[0][qq][1] << "  qqqqqqq   " << finalestimate[0][qq][0] << " " << finalestimate[0][qq][1] << "      " <<finalestimate2[0][qq][0] << " " << finalestimate2[0][qq][1]<< endl;
			//}
            //cout << "happy" <<endl;
			//for(int qq = 0;qq<6;++qq){
			//	cout << chCoef[1][qq][0] << " " << chCoef[1][qq][1]<< "         2bb       " << chCoef2[1][qq][0] << " " << chCoef2[1][qq][1] << "  qqqqqqq   " << finalestimate[1][qq][0] << " " << finalestimate[1][qq][1] << "      " <<finalestimate2[1][qq][0] << " " << finalestimate2[1][qq][1]<< endl;
			//}
			fourwayMLDT(pow(stdDev, 2), chCoef, chCoef2, rx, app, appLlr, finalestimate,finalestimate2, scma_matrix, coef_idx);
			Detecter(data, appLlr, error);
			ber[i] = error / ((long double)NUM_USER * (SCMA_SOURCE / SCMA_USER_SOURCE) * block * (BLOCK_LEN - DIFF_ENC) * Qm); // excluding reference symbol
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

