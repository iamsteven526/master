#define _USE_MATH_DEFINES
#include <iostream>
#include <random>
#include "parameters.h"
#include <cstring>
using namespace std;
#pragma warning(disable:4996)
namespace
{
	random_device seed;
	mt19937 generator(seed());
	uniform_real_distribution<double> uniform(0, 1);
	normal_distribution<double> normal(0, 1);
}

void Clustering(double *****postRx, double variance, double *****estimate, double *****H, int ***cluster_num, int *packet_num, double& mse, bool** estimate_packet_time, double** cluster_sample, double** centroid, int** group, int* groupSize, double *variation, int** packet_time, double **softAssign)
{
	
	if (SLOTED)
	{
		double total_mse = 0;
		int total_num = 0;
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			for (int t = 0; t < packet_num[nuser]; t++)
			{
				for (int i = 0; i < FFT_POINT; i++)
				{
					estimate[t][nuser][0][i][0] = H[t][nuser][0][0][i] + sqrt(variance / (2 * (FFT_SEGMENT + DIFF_ENC))) * normal(generator);
					estimate[t][nuser][0][i][1] = H[t][nuser][0][1][i] + sqrt(variance / (2 * (FFT_SEGMENT + DIFF_ENC))) * normal(generator);
					//estimate[t][nuser][0][0][i] = H[t][nuser][0][0][i] + sqrt(variance / c[nuser]) * normal(generator);
					//estimate[t][nuser][0][1][i] = H[t][nuser][0][1][i] + sqrt(variance / c[nuser]) * normal(generator);
					total_mse += abs(sqrt(variance / (2 * (FFT_SEGMENT + DIFF_ENC))) * normal(generator));
					total_num++;
				}
			}
		}
		if(total_num!=0)
			mse += total_mse / total_num;
	}
	else
	{
		if (CH_TYPE == 1)
		{
			//int cluster_len = 60;
			int cluster_len = (packet_time[1][0] - packet_time[0][0]) / (FFT_POINT + CP_LEN) - PREABLE_LEN / (FFT_POINT + CP_LEN) - 1;
			int comple = FFT_SEGMENT - cluster_len - PREABLE_LEN / (FFT_POINT + CP_LEN);			
			if (CE_METHOD == 1 && cluster_len > 0 && cluster_len < 61)
			{
				for (int nuser = 0; nuser < NUM_USER; nuser++)
				{
					for (int t = 0; t < packet_num[nuser]; t++)
					{
						if (Time_Estimate && estimate_packet_time[nuser][t] == 0)
							continue;
						
						for (int i = 0; i < FFT_POINT; i++)
						{
							bool flag = 1;
							int count = 0;
							centroid[0][0] = estimate[t][nuser][0][i][0];
							centroid[0][1] = estimate[t][nuser][0][i][1];
							centroid[1][0] = -estimate[t][nuser][0][i][0];
							centroid[1][1] = -estimate[t][nuser][0][i][1];


							for (int j = 0; j < cluster_len; j++)
							{
								cluster_sample[j][0] = postRx[t][nuser][j + nuser * comple][0][i];
								cluster_sample[j][1] = postRx[t][nuser][j + nuser * comple][1][i];
								cluster_sample[cluster_len + j][0] = -cluster_sample[j][0];
								cluster_sample[cluster_len + j][1] = -cluster_sample[j][1];

								// weight of preamble estimation
								cluster_sample[2 * cluster_len][0] = estimate[t][nuser][0][i][0];
								cluster_sample[2 * cluster_len][1] = estimate[t][nuser][0][i][1];
								cluster_sample[2 * cluster_len + 1][0] = -estimate[t][nuser][0][i][0];
								cluster_sample[2 * cluster_len + 1][1] = -estimate[t][nuser][0][i][1];

								cluster_sample[2 * cluster_len + 2][0] = estimate[t][nuser][0][i][0];
								cluster_sample[2 * cluster_len + 2][1] = estimate[t][nuser][0][i][1];
								cluster_sample[2 * cluster_len + 3][0] = -estimate[t][nuser][0][i][0];
								cluster_sample[2 * cluster_len + 3][1] = -estimate[t][nuser][0][i][1];
							}

							

							/*if (EM_GMM&& flag==1) EMClustering(cluster_sample, centroid, softAssign, variance, 2 * cluster_len + 4, 2);
							else
							{
								Grouping(cluster_sample, centroid, group, groupSize, 1, 2, 2 * cluster_len + 4);
								CentroidRenewal(cluster_sample, centroid, group, groupSize, variation, 2);
							}

							if (centroid[0][0] < 0.001 && centroid[0][1] < 0.001)
							{
								flag = 0;
								continue;
							}*/

							Grouping(cluster_sample, centroid, group, groupSize, 1, 2, 2 * cluster_len + 4);
							CentroidRenewal(cluster_sample, centroid, group, groupSize, variation, 2);
                            
							double temp = pow(H[t][nuser][0][0][i] - centroid[0][0], 2) + pow(H[t][nuser][0][1][i] - centroid[0][1], 2) - pow(H[t][nuser][0][0][i] + centroid[0][0], 2) - pow(H[t][nuser][0][1][i] + centroid[0][1], 2);

							estimate[t][nuser][0][i][0] = temp < 0 ? centroid[0][0] : -centroid[0][0];
							estimate[t][nuser][0][i][1] = temp < 0 ? centroid[0][1] : -centroid[0][1];
							

							for (int j = 1; j < FFT_SEGMENT; j++)
							{
								estimate[t][nuser][j][i][0] = estimate[t][nuser][0][i][0];
								estimate[t][nuser][j][i][1] = estimate[t][nuser][0][i][1];
							}
						}
					}
				}
			}


			double total_mse = 0;
			int total_num = 0;
			for (int nuser = 0; nuser < NUM_USER; nuser++)
			{
				for (int t = 0; t < packet_num[nuser]; t++)
				{
					if (Time_Estimate && estimate_packet_time[nuser][t] == 0)
						continue;

					for (int i = 0; i < FFT_POINT; i++)
					{
						double real = H[0][nuser][0][0][i] - estimate[t][nuser][0][i][0];
						double imag = H[0][nuser][0][1][i] - estimate[t][nuser][0][i][1];

						total_mse += 0.5 * (pow(real, 2) + pow(imag, 2));
						total_num++;
					}
				}
			}
			if (total_num != 0)
				mse += total_mse / total_num;
		}
		else if (CH_TYPE == 2)
		{
			//---- set preamble estimation gain
			for (int nuser = 0; nuser < NUM_USER; nuser++)
			{
				for (int t = 0; t < packet_num[nuser]; t++)
				{
					if (Time_Estimate && estimate_packet_time[nuser][t] == 0)
						continue;

					for (int i = 0; i < FFT_POINT; i++)
					{
						for (int j = 1; j < FFT_SEGMENT; j++)
						{
							estimate[t][nuser][j][i][0] = estimate[t][nuser][0][i][0];
							estimate[t][nuser][j][i][1] = estimate[t][nuser][0][i][1];
						}
					}
				}
			}

			int cluster_len = (packet_time[1][0] - packet_time[0][0]) / (FFT_POINT + CP_LEN) - 1;
			//int cluster_len = 60;
			//cout << cluster_len << endl;
			
			if (CE_METHOD == 1  && cluster_len < 61)
			{
				if (cluster_len > WINDOW_SIZE)
				{
					//-------------user 1
					for (int w = 0; w < cluster_len - WINDOW_SIZE; w++)
					{
						for (int i = 0; i < FFT_POINT; i++)
						{
							
							centroid[0][0] = estimate[0][0][w + WINDOW_SIZE - 7][i][0];
							centroid[0][1] = estimate[0][0][w + WINDOW_SIZE - 7][i][1];
							centroid[1][0] = -estimate[0][0][w + WINDOW_SIZE - 7][i][0];
							centroid[1][1] = -estimate[0][0][w + WINDOW_SIZE - 7][i][1];


							for (int j = 0; j < WINDOW_SIZE; j++)
							{
								cluster_sample[j][0] = postRx[0][0][w + j][0][i];
								cluster_sample[j][1] = postRx[0][0][w + j][1][i];
								cluster_sample[WINDOW_SIZE + j][0] = -cluster_sample[j][0];
								cluster_sample[WINDOW_SIZE + j][1] = -cluster_sample[j][1];
							}

							Grouping(cluster_sample, centroid, group, groupSize, 1, 2, 2 * WINDOW_SIZE);
							CentroidRenewal(cluster_sample, centroid, group, groupSize, variation, 2);

							double temp = pow(H[0][0][w + WINDOW_SIZE - 6][0][i] - centroid[0][0], 2) + pow(H[0][0][w + WINDOW_SIZE - 6][1][i] - centroid[0][1], 2) - pow(H[0][0][w + WINDOW_SIZE - 6][0][i] + centroid[0][0], 2) - pow(H[0][0][w + WINDOW_SIZE - 6][1][i] + centroid[0][1], 2);
							
							//estimate[0][0][w + WINDOW_SIZE - 6][i][0] = temp < 0 ? centroid[0][0] : -centroid[0][0];
							//estimate[0][0][w + WINDOW_SIZE - 6][i][1] = temp < 0 ? centroid[0][1] : -centroid[0][1];

							estimate[0][0][w + WINDOW_SIZE - 6][i][0] = H[0][0][w + WINDOW_SIZE - 6][0][i];
							estimate[0][0][w + WINDOW_SIZE - 6][i][1] = H[0][0][w + WINDOW_SIZE - 6][1][i];
						}
					}

					for (int j = cluster_len - 6; j < FFT_SEGMENT; j++)
					{
						for (int i = 0; i < FFT_POINT; i++)
						{
							estimate[0][0][j][i][0] = estimate[0][0][cluster_len - 7][i][0];
							estimate[0][0][j][i][1] = estimate[0][0][cluster_len - 7][i][1];
						}
					}

					for (int j = 0; j < 2 * FFT_SEGMENT; j++)
					{
						cluster_sample[j][0] = 0;
						cluster_sample[j][1] = 0;
					}
					//-------------user 2
					for (int w = 0; w < cluster_len - WINDOW_SIZE; w++)
					{
						for (int i = 0; i < FFT_POINT; i++)
						{
							int count = 0;
							centroid[0][0] = H[0][1][FFT_SEGMENT - (WINDOW_SIZE + w - 5) ][0][i];
							centroid[0][1] = H[0][1][FFT_SEGMENT - (WINDOW_SIZE + w - 5)][1][i];
							centroid[1][0] = - centroid[0][0];
							centroid[1][1] = - centroid[0][1];


							for (int j = 0; j < WINDOW_SIZE; j++)
							{
								cluster_sample[j][0] = postRx[0][1][FFT_SEGMENT - j - 1 - w][0][i];
								cluster_sample[j][1] = postRx[0][1][FFT_SEGMENT - j - 1 - w][1][i];
								cluster_sample[WINDOW_SIZE + j][0] = -cluster_sample[j][0];
								cluster_sample[WINDOW_SIZE + j][1] = -cluster_sample[j][1];
							}

							Grouping(cluster_sample, centroid, group, groupSize, 1, 2, 2 * WINDOW_SIZE);
							CentroidRenewal(cluster_sample, centroid, group, groupSize, variation, 2);

							double temp = pow(H[0][1][FFT_SEGMENT - (WINDOW_SIZE + w - 5)][0][i] - centroid[0][0], 2) + pow(H[0][1][FFT_SEGMENT - (WINDOW_SIZE + w - 5)][1][i] - centroid[0][1], 2) - pow(H[0][1][FFT_SEGMENT - (WINDOW_SIZE + w - 5)][0][i] + centroid[0][0], 2) - pow(H[0][1][FFT_SEGMENT - (WINDOW_SIZE + w - 5)][1][i] + centroid[0][1], 2);
							
							//estimate[0][1][FFT_SEGMENT - (WINDOW_SIZE + w - 5)][i][0] = temp < 0 ? centroid[0][0] : -centroid[0][0];
							//estimate[0][1][FFT_SEGMENT - (WINDOW_SIZE + w - 5)][i][1] = temp < 0 ? centroid[0][1] : -centroid[0][1];
							estimate[0][1][FFT_SEGMENT - (WINDOW_SIZE + w - 5)][i][0] = H[0][1][FFT_SEGMENT - (WINDOW_SIZE + w - 5)][0][i];
							estimate[0][1][FFT_SEGMENT - (WINDOW_SIZE + w - 5)][i][1] = H[0][1][FFT_SEGMENT - (WINDOW_SIZE + w - 5)][1][i];
						}
					}
					//system("pause");

					//cout << "****"<<estimate[0][1][FFT_SEGMENT - WINDOW_SIZE - 1][0][0] << "****" << endl;
					for (int j = 0; j < WINDOW_SIZE / 2; j++)
					{
						for (int i = 0; i < FFT_POINT; i++)
						{	
							estimate[0][1][FFT_SEGMENT - j - 1][i][0] = estimate[0][1][FFT_SEGMENT - WINDOW_SIZE/2 - 1][i][0];
							estimate[0][1][FFT_SEGMENT - j - 1][i][1] = estimate[0][1][FFT_SEGMENT - WINDOW_SIZE/2 - 1][i][1];
						}
					}

					for (int j = 0; j < 2 * FFT_SEGMENT; j++)
					{
						cluster_sample[j][0] = 0;
						cluster_sample[j][1] = 0;
					}
					
				}
			}
		}
	}

	/*for (int i = 0; i < FFT_SEGMENT; i++)
	{
		cout << H[0][0][i][0][1] << " ";
	}
	cout << endl;

	for (int i = 0; i < FFT_SEGMENT; i++)
	{
		cout << estimate[0][0][i][1][0] << " ";
	}
	//system("pause");

	cout << endl;
	for (int i = 0; i < FFT_SEGMENT; i++)
	{
		cout << H[0][1][i][0][0] << " ";
	}
	cout << endl;

	for (int i = 0; i < FFT_SEGMENT; i++)
	{
		cout << estimate[0][1][i][0][0] << " ";
	}
	system("pause");*/
}

void Grouping(double **rx, double **centroid, int **group, int *groupSize, int CLUSTER_USER, int CLUSTER_GROUP, int CLUSTER_LEN)
{
	memset(groupSize, 0, sizeof(int) * CLUSTER_GROUP);
	vector<double> dist(CLUSTER_GROUP);

	for (int i = 0; i < CLUSTER_LEN; i++)
	{
		for (int j = 0; j < CLUSTER_GROUP; j++)
		{
			dist[j] = EuclideanDistance(rx[i], centroid[j]);
		}

		int min = 0;
		for (int j = 1; j < CLUSTER_GROUP; j++)
		{
			if (dist[j] < dist[min]) min = j;
		}
		group[min][groupSize[min]++] = i;
	}

}

double EuclideanDistance(double* x, double* y)
{
	return sqrt(pow(x[0] - y[0], 2) + pow(x[1] - y[1], 2));
}

void CentroidRenewal(double **rx, double **centroid, int **group, int *groupSize, double *variation, int CLUSTER_GROUP)
{
	memset(variation, 0, sizeof(double) * CLUSTER_GROUP);
	for (int i = 0; i < CLUSTER_GROUP; i++)
	{
		for (int j = 0; j < 2; j++) // real and imaginary
		{
			double temp = 0;
			for (int k = 0; k < groupSize[i]; k++)
			{
				temp += rx[group[i][k]][j];
			}
			temp /= groupSize[i]; //�s��centroid
			//cout << groupSize[i] << " ";
			variation[i] += pow(centroid[i][j] - temp, 2);	//�s�¤������t�O
			centroid[i][j] = temp;
			//cout << centroid[i][j] << " ";
		}
		variation[i] = sqrt(variation[i]);
	}
//	cout << endl;
}

bool ConditionCheck(double** centroid, int* groupSize, int CLUSTER_GROUP)
{
	for (int i = 0; i < CLUSTER_GROUP; i++)//�ˬdGROUP_SIZE��centroid�����I
	{
		if (groupSize[i] == 0) return true; // empty cluster
	}
	double average = 0;
	for (int i = 0; i < 2; i++)
	{
		double temp = 0;
		for (int j = 0; j < CLUSTER_GROUP; j++)
		{
			temp += centroid[j][i];
		}
		average += pow(temp / CLUSTER_GROUP, 2);
	}
	if (sqrt(average) > 0.05) // ad-hoc
	{
		return true;
	}
	return false;
}

void EMClustering(double** rx, double** centroid, double** softAssign, double variance, int CLUSTER_LEN, int CLUSTER_GROUP)
{
	int it = 0;
	while (it++ < 1000)
	{
		for (int i = 0; i < CLUSTER_LEN; i++)
		{
			double sum = 0;
			for (int k = 0; k < CLUSTER_GROUP; k++)
			{
				softAssign[i][k] = exp(-pow(rx[i][0] - centroid[k][0], 2) / (2. * variance)) * exp(-pow(rx[i][1] - centroid[k][1], 2) / (2. * variance));
				if (softAssign[i][k] < NUMERIC_LIMIT) softAssign[i][k] = NUMERIC_LIMIT;
				sum += softAssign[i][k];
			}
			for (int k = 0; k < CLUSTER_GROUP; k++)
			{
				softAssign[i][k] /= sum;
				if (softAssign[i][k] < NUMERIC_LIMIT) softAssign[i][k] = NUMERIC_LIMIT;
			}
		}
		double variation = 0;
		for (int k = 0; k < CLUSTER_GROUP; k++)
		{
			double eMean[2] = { 0 }, nFactor = 0;
			for (int i = 0; i < CLUSTER_LEN; i++)
			{
				eMean[0] += softAssign[i][k] * rx[i][0];
				eMean[1] += softAssign[i][k] * rx[i][1];
				nFactor += softAssign[i][k];
			}
			eMean[0] /= nFactor;
			eMean[1] /= nFactor;
			variation += pow(centroid[k][0] - eMean[0], 2) + pow(centroid[k][1] - eMean[1], 2);
			centroid[k][0] = eMean[0]; centroid[k][1] = eMean[1];
		}
		if (variation < 1e-10) return;
	}
}
