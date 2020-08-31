#define _USE_MATH_DEFINES
#include <iostream>
#include <random>
#include "parameters.h"
using namespace std;
#pragma warning(disable:4996)
namespace
{
	random_device seed;
	mt19937 generator(seed());
	uniform_real_distribution<double> uniform(0, 1);
}

void Clustering(double ***postRx, double **centroid, int **group, int *groupSize, double *distList, double *variation, double **softAssign, double variance, double ****estimate)
{
	
	/*for (int sch = 0; sch < FFT_POINT; sch++) // estimation is performed for each sub-channel
	{
		//---- fliping
		for (int i = 0; i < FFT_SEGMENT + DIFF_ENC; i++)
		{
			postRx[i + FFT_SEGMENT + DIFF_ENC][0][sch] = -postRx[i][0][sch];
			postRx[i + FFT_SEGMENT + DIFF_ENC][1][sch] = -postRx[i][1][sch];
		}

		bool reset = true; // resultant centroids will become the initial centroids in the next windowing
		for (int offset = 0; offset <= (FFT_SEGMENT + DIFF_ENC - WINDOW_SIZE)*SLIDING; offset += 1) // sliding window
		{
			//bool reset = true; // select a new initial centroids in each windowing
			int count = 0;
			while (true)
			{
				if (reset)
				{
					InitialSeeding(postRx, centroid, group, groupSize, distList, sch, offset);
					reset = false;
				}
				Grouping(postRx, centroid, group, groupSize, sch, offset);
				CentroidRenewal(postRx, centroid, group, groupSize, variation, sch);
				bool convergence = true;
				for (int i = 0; i < GROUP_SIZE; i++)
				{
					if (variation[i] > 1e-5)
					{
						convergence = false;
						break;
					}
				}
				if (convergence)
				{
					reset = ConditionCheck(centroid, groupSize, count);
					if (!reset)
					{
						if (EM_GMM) EMClustering(postRx, centroid, softAssign, variance, sch, offset);
						CoefEstimation(centroid, estimate, sch, offset);
						break;
					}
				}
			}
		}
	}*/
}

void InitialSeeding(double ***postRx, double **centroid, int **group, int *groupSize, double *distList, int sch, int offset)
{
	if (INI_METHOD == 0) // LBG
	{
		double delta = 0.001;
		double phi = 2 * M_PI * uniform(generator); // splitting phase
		double dist[GROUP_SIZE];
		centroid[0][0] = +delta*cos(phi); centroid[0][1] = +delta*sin(phi);
		centroid[1][0] = -delta*cos(phi); centroid[1][1] = -delta*sin(phi);
		for (int i = 1; i < NUM_USER; i++)
		{
			memset(groupSize, 0, sizeof(int) * GROUP_SIZE);
			for (int j = offset; j < (SLIDING ? offset + WINDOW_SIZE : FFT_SEGMENT + DIFF_ENC); j++) // grouping
			{
				for (int k = 0; k < (1 << i); k++)
				{
					dist[k] = sqrt(pow(postRx[j][0][sch] - centroid[k][0], 2) + pow(postRx[j][1][sch] - centroid[k][1], 2));
				}
				int min = 0;
				for (int k = 1; k < (1 << i); k++)
				{
					if (dist[k] < dist[min]) min = k;
				}
				group[min][groupSize[min]++] = j;
			}
			for (int j = 0; j < (1 << i); j++) // splitting
			{
				centroid[2 * j][0] = centroid[2 * j][1] = 0;
				for (int k = 0; k < groupSize[j]; k++)
				{
					centroid[2 * j][0] += postRx[group[j][k]][0][sch];
					centroid[2 * j][1] += postRx[group[j][k]][1][sch];
				}
				centroid[2 * j + 0][0] /= groupSize[j]; centroid[2 * j + 1][0] = centroid[2 * j + 0][0];
				centroid[2 * j + 0][1] /= groupSize[j]; centroid[2 * j + 1][1] = centroid[2 * j + 0][1];
				centroid[2 * j + 0][0] += delta*cos(phi); centroid[2 * j + 0][1] += delta*sin(phi);
				centroid[2 * j + 1][0] -= delta*cos(phi); centroid[2 * j + 1][1] -= delta*sin(phi);
			}
		}
	}
	else if (INI_METHOD == 1 || INI_METHOD == 2) // k-means++
	{
		int iniPick = offset + int(rand()% (2*(FFT_SEGMENT + DIFF_ENC)));
		
			centroid[0][0] = postRx[iniPick][0][sch]; centroid[0][1] = postRx[iniPick][1][sch];
		for (int i = 1; i < GROUP_SIZE; i++)
		{
			double sum = 0;
			for (int j = offset; j < (SLIDING ? offset + WINDOW_SIZE : FFT_SEGMENT + DIFF_ENC); j++) // distance from the closest centroid
			{
				distList[j] = sqrt(pow(postRx[j][0][sch] - centroid[0][0], 2) + pow(postRx[j][1][sch] - centroid[0][1], 2));
				for (int k = 1; k < i; k++)
				{
					double temp = sqrt(pow(postRx[j][0][sch] - centroid[k][0], 2) + pow(postRx[j][1][sch] - centroid[k][1], 2));
					if (temp < distList[j]) distList[j] = temp;
				}
				distList[j] = pow(distList[j], 2);
				sum += distList[j];
			}
			if (INI_METHOD == 2)
			{
				double list[2][SLIDING ? WINDOW_SIZE : FFT_SEGMENT + DIFF_ENC];
				memcpy(list[1], &distList[offset], sizeof(double)*(SLIDING ? WINDOW_SIZE : FFT_SEGMENT + DIFF_ENC));
				for (int j = 0; j < (SLIDING ? WINDOW_SIZE : FFT_SEGMENT + DIFF_ENC); j++)
				{
					list[0][j] = offset + j;
				}
				for (int j = 0; j < (SLIDING ? WINDOW_SIZE : FFT_SEGMENT + DIFF_ENC) - 1; j++) // sorting
				{
					for (int k = (SLIDING ? WINDOW_SIZE : FFT_SEGMENT + DIFF_ENC) - 1; k > j; k--)
					{
						if (list[1][k] < list[1][k - 1])
						{
							double temp = list[0][k];
							list[0][k] = list[0][k - 1];
							list[0][k - 1] = temp;
							temp = list[1][k];
							list[1][k] = list[1][k - 1];
							list[1][k - 1] = temp;
						}
					}
				}
				for (int j = 0; j < (SLIDING ? WINDOW_SIZE : FFT_SEGMENT + DIFF_ENC) - ((SLIDING ? WINDOW_SIZE : FFT_SEGMENT + DIFF_ENC) / GROUP_SIZE); j++)
				{
					sum -= distList[(int)list[0][j]];
					distList[(int)list[0][j]] = 0;
				}
			}
			double ranNum = uniform(generator)*sum;
			for (int j = offset; j < (SLIDING ? offset + WINDOW_SIZE : FFT_SEGMENT + DIFF_ENC); j++) // select with weighted probability
			{
				ranNum -= distList[j];
				if (ranNum <= 0)
				{
					centroid[i][0] = postRx[j][0][sch];
					centroid[i][1] = postRx[j][1][sch];
					break;
				}
			}
		}
	}
	else
	{
		printf("\nPARAMETER SETTING IS WRONG\n");
		system("pause");
	}
}

void Grouping(double ***postRx, double **centroid, int **group, int *groupSize, int sch, int offset)
{
	memset(groupSize, 0, sizeof(int) * GROUP_SIZE);
	double dist[GROUP_SIZE];

	if (SLIDING)
	{
		int i = 0;
		for (int k = 0 ; k <  2 * WINDOW_SIZE ; k++)
		{
			if (k < WINDOW_SIZE)
				i = k + offset;
			else
				i = k + offset + FFT_SEGMENT + DIFF_ENC - WINDOW_SIZE;
			
			
			for (int j = 0; j < GROUP_SIZE; j++)
			{
				dist[j] = sqrt(pow(postRx[i][0][sch] - centroid[j][0], 2) + pow(postRx[i][1][sch] - centroid[j][1], 2));
			}
			int min = 0;
			for (int j = 1; j < GROUP_SIZE; j++)
			{
				if (dist[j] < dist[min]) min = j;
			}
			group[min][groupSize[min]++] = i;
		}
	}
	else
	{
		for (int i = offset; i < 2 * (FFT_SEGMENT + DIFF_ENC); i++)
		{
			for (int j = 0; j < GROUP_SIZE; j++)
			{
				dist[j] = sqrt(pow(postRx[i][0][sch] - centroid[j][0], 2) + pow(postRx[i][1][sch] - centroid[j][1], 2));
			}
			int min = 0;
			for (int j = 1; j < GROUP_SIZE; j++)
			{
				if (dist[j] < dist[min]) min = j;
			}
			group[min][groupSize[min]++] = i;
		}
	}
	
}

void CentroidRenewal(double ***postRx, double **centroid, int **group, int *groupSize, double *variation, int sch)
{
	memset(variation, 0, sizeof(double) * GROUP_SIZE);
	for (int i = 0; i < GROUP_SIZE; i++)
	{
		for (int j = 0; j < 2; j++) // real and imaginary
		{
			double temp = 0;
			for (int k = 0; k < groupSize[i]; k++)
			{
				temp += postRx[group[i][k]][j][sch];
			}
			temp /= groupSize[i];
			variation[i] += pow(centroid[i][j] - temp, 2);
			centroid[i][j] = temp;
		}
		variation[i] = sqrt(variation[i]);
	}
}

void CoefEstimation(double **centroid, double ****estimate, int sch, int offset)
{
	if (NUM_USER == 1)
	{
		estimate[0][(offset + (WINDOW_SIZE - 1) / 2)*SLIDING][sch][0] = (centroid[0][0] - centroid[1][0])*0.5;
		estimate[0][(offset + (WINDOW_SIZE - 1) / 2)*SLIDING][sch][1] = (centroid[0][1] - centroid[1][1])*0.5;
	}
	else if (NUM_USER == 2)
	{
		double dev = pow(centroid[0][0] + centroid[1][0], 2) + pow(centroid[0][1] + centroid[1][1], 2);
		int index = 1, count = 0;
		for (int i = 2; i < GROUP_SIZE; i++)
		{
			double temp = pow(centroid[0][0] + centroid[i][0], 2) + pow(centroid[0][1] + centroid[i][1], 2);
			if (temp < dev)
			{
				dev = temp;
				index = i;
			}
		}
		estimate[0][(offset + (WINDOW_SIZE - 1) / 2)*SLIDING][sch][0] = estimate[0][(offset + (WINDOW_SIZE - 1) / 2)*SLIDING][sch][1] = 0;
		estimate[1][(offset + (WINDOW_SIZE - 1) / 2)*SLIDING][sch][0] = estimate[1][(offset + (WINDOW_SIZE - 1) / 2)*SLIDING][sch][1] = 0;
		for (int i = 1; i < GROUP_SIZE; i++)
		{
			if (i != index)
			{
				estimate[count][(offset + (WINDOW_SIZE - 1) / 2)*SLIDING][sch][0] += 0.25*(centroid[0][0] + centroid[i][0]);
				estimate[count][(offset + (WINDOW_SIZE - 1) / 2)*SLIDING][sch][1] += 0.25*(centroid[0][1] + centroid[i][1]);
				estimate[(count + 1) % NUM_USER][(offset + (WINDOW_SIZE - 1) / 2)*SLIDING][sch][0] -= 0.25*(centroid[index][0] + centroid[i][0]);
				estimate[(count + 1) % NUM_USER][(offset + (WINDOW_SIZE - 1) / 2)*SLIDING][sch][1] -= 0.25*(centroid[index][1] + centroid[i][1]);
				count++;
			}
		}
	}
	else if (NUM_USER == 3)
	{
		int pair[GROUP_SIZE >> 1][2], flag[GROUP_SIZE] = { 0 }, count = 0;
		double tCentroid[GROUP_SIZE >> 1][2], midpoint[GROUP_SIZE][2];
		double **sum = new double *[GROUP_SIZE - 1];
		for (int i = 0; i < GROUP_SIZE - 1; i++)
		{
			sum[i] = new double[GROUP_SIZE - i - 1];
			for (int j = 0; j < GROUP_SIZE - i - 1; j++)
			{
				sum[i][j] = pow(centroid[i][0] + centroid[i + j + 1][0], 2) + pow(centroid[i][1] + centroid[i + j + 1][1], 2);
			}
		}
		while (count < (GROUP_SIZE >> 1)) // pairing
		{
			double min = 1e10;
			for (int i = 0; i < GROUP_SIZE - 1; i++)
			{
				if (flag[i] == 0)
				{
					for (int j = 0; j < GROUP_SIZE - i - 1; j++)
					{
						if (sum[i][j] < min && flag[i + j + 1] == 0)
						{
							min = sum[i][j];
							pair[count][0] = i;
							pair[count][1] = i + j + 1;
						}
					}
				}
			}
			flag[pair[count][0]] = flag[pair[count][1]] = 1;
			count++;
		}
		for (int i = 1; i < (GROUP_SIZE >> 1); i++)
		{
			tCentroid[0][0] = centroid[pair[0][0]][0]; tCentroid[0][1] = centroid[pair[0][0]][1];
			tCentroid[1][0] = centroid[pair[0][1]][0]; tCentroid[1][1] = centroid[pair[0][1]][1];
			tCentroid[2][0] = centroid[pair[i][0]][0]; tCentroid[2][1] = centroid[pair[i][0]][1];
			tCentroid[3][0] = centroid[pair[i][1]][0]; tCentroid[3][1] = centroid[pair[i][1]][1];
			double dev = pow(tCentroid[0][0] + tCentroid[1][0], 2) + pow(tCentroid[0][1] + tCentroid[1][1], 2);
			int index = 1, count = 0;
			for (int j = 2; j < (GROUP_SIZE >> 1); j++)
			{
				double temp = pow(tCentroid[0][0] + tCentroid[j][0], 2) + pow(tCentroid[0][1] + tCentroid[j][1], 2);
				if (temp < dev)
				{
					dev = temp;
					index = j;
				}
			}
			for (int j = 1; j < (GROUP_SIZE >> 1); j++)
			{
				if (j != index)
				{
					midpoint[count][0] = 0.5*(tCentroid[0][0] + tCentroid[j][0]);
					midpoint[count][1] = 0.5*(tCentroid[0][1] + tCentroid[j][1]);
					midpoint[count + 1][0] = 0.5*(tCentroid[index][0] + tCentroid[j][0]);
					midpoint[count + 1][1] = 0.5*(tCentroid[index][1] + tCentroid[j][1]);
					count += 2;
				}
			}
			tCentroid[0][0] = centroid[pair[((i + 0) % ((GROUP_SIZE >> 1) - 1)) + 1][0]][0]; tCentroid[0][1] = centroid[pair[((i + 0) % ((GROUP_SIZE >> 1) - 1)) + 1][0]][1];
			tCentroid[1][0] = centroid[pair[((i + 0) % ((GROUP_SIZE >> 1) - 1)) + 1][1]][0]; tCentroid[1][1] = centroid[pair[((i + 0) % ((GROUP_SIZE >> 1) - 1)) + 1][1]][1];
			tCentroid[2][0] = centroid[pair[((i + 1) % ((GROUP_SIZE >> 1) - 1)) + 1][0]][0]; tCentroid[2][1] = centroid[pair[((i + 1) % ((GROUP_SIZE >> 1) - 1)) + 1][0]][1];
			tCentroid[3][0] = centroid[pair[((i + 1) % ((GROUP_SIZE >> 1) - 1)) + 1][1]][0]; tCentroid[3][1] = centroid[pair[((i + 1) % ((GROUP_SIZE >> 1) - 1)) + 1][1]][1];
			dev = pow(tCentroid[0][0] + tCentroid[1][0], 2) + pow(tCentroid[0][1] + tCentroid[1][1], 2);
			index = 1;
			for (int j = 2; j < (GROUP_SIZE >> 1); j++)
			{
				double temp = pow(tCentroid[0][0] + tCentroid[j][0], 2) + pow(tCentroid[0][1] + tCentroid[j][1], 2);
				if (temp < dev)
				{
					dev = temp;
					index = j;
				}
			}
			for (int j = 1; j < (GROUP_SIZE >> 1); j++)
			{
				if (j != index)
				{
					midpoint[count][0] = 0.5*(tCentroid[0][0] + tCentroid[j][0]);
					midpoint[count][1] = 0.5*(tCentroid[0][1] + tCentroid[j][1]);
					midpoint[count + 1][0] = 0.5*(tCentroid[index][0] + tCentroid[j][0]);
					midpoint[count + 1][1] = 0.5*(tCentroid[index][1] + tCentroid[j][1]);
					count += 2;
				}
			}
			for (int j = 0; j < 2; j++)
			{
				double temp = midpoint[j*(GROUP_SIZE >> 1) + 2][0];
				midpoint[j*(GROUP_SIZE >> 1) + 2][0] = midpoint[j*(GROUP_SIZE >> 1) + 3][0];
				midpoint[j*(GROUP_SIZE >> 1) + 3][0] = temp;
				temp = midpoint[j*(GROUP_SIZE >> 1) + 2][1];
				midpoint[j*(GROUP_SIZE >> 1) + 2][1] = midpoint[j*(GROUP_SIZE >> 1) + 3][1];
				midpoint[j*(GROUP_SIZE >> 1) + 3][1] = temp;
			}
			double min = 1e10;
			int target[2];
			for (int j = 0; j < (GROUP_SIZE >> 1); j++) // find the closest pair
			{
				for (int k = 0; k < (GROUP_SIZE >> 1); k++)
				{
					double temp = pow(midpoint[j][0] - midpoint[(GROUP_SIZE >> 1) + k][0], 2) + pow(midpoint[j][1] - midpoint[(GROUP_SIZE >> 1) + k][1], 2);
					if (temp < min)
					{
						min = temp;
						target[0] = j;
						target[1] = k;
					}
				}
			}
			estimate[i - 1][(offset + (WINDOW_SIZE - 1) / 2)*SLIDING][sch][0] = (midpoint[target[0]][0] - midpoint[(target[0] + 2) % (GROUP_SIZE >> 1)][0]) / 4.;
			estimate[i - 1][(offset + (WINDOW_SIZE - 1) / 2)*SLIDING][sch][1] = (midpoint[target[0]][1] - midpoint[(target[0] + 2) % (GROUP_SIZE >> 1)][1]) / 4.;
			estimate[i - 1][(offset + (WINDOW_SIZE - 1) / 2)*SLIDING][sch][0] += (midpoint[(GROUP_SIZE >> 1) + target[1]][0] - midpoint[(GROUP_SIZE >> 1) + ((target[1] + 2) % (GROUP_SIZE >> 1))][0]) / 4.;
			estimate[i - 1][(offset + (WINDOW_SIZE - 1) / 2)*SLIDING][sch][1] += (midpoint[(GROUP_SIZE >> 1) + target[1]][1] - midpoint[(GROUP_SIZE >> 1) + ((target[1] + 2) % (GROUP_SIZE >> 1))][1]) / 4.;
		}
		for (int i = 0; i < GROUP_SIZE - 1; i++)
		{
			delete[] sum[i];
		}
		delete[] sum;
	}
}

void UserIdentification(double ****H, double ****estimate)
{
	if (CH_TYPE == 2 && !SLIDING) // average the channel coefficients varing with time for comparison 
	{
		for (int i = 0; i < NUM_USER; i++)
		{
			for (int j = 1; j < FFT_SEGMENT + DIFF_ENC; j++)
			{
				for (int k = 0; k < FFT_POINT; k++)
				{
					H[i][0][0][k] += H[i][j][0][k];
					H[i][0][1][k] += H[i][j][1][k];
				}
			}
			for (int k = 0; k < FFT_POINT; k++)
			{
				H[i][0][0][k] /= (double)(FFT_SEGMENT + DIFF_ENC);
				H[i][0][1][k] /= (double)(FFT_SEGMENT + DIFF_ENC);
			}
		}
	}

	for (int j = ((WINDOW_SIZE - 1) / 2)*SLIDING; j <= ((WINDOW_SIZE - 1) / 2 + (FFT_SEGMENT + DIFF_ENC - WINDOW_SIZE))*SLIDING; j++)
	{
		for (int sch = 0; sch < FFT_POINT; sch++)
		{
			if (NUM_USER == 2)
			{
				double comparison[2] = { 0 };
				for (int i = 0; i < NUM_USER; i++)
				{
					double temp1 = pow(H[i][j][0][sch] - estimate[i][j][sch][0], 2) + pow(H[i][j][1][sch] - estimate[i][j][sch][1], 2);
					double temp2 = pow(H[i][j][0][sch] + estimate[i][j][sch][0], 2) + pow(H[i][j][1][sch] + estimate[i][j][sch][1], 2);
					temp1 < temp2 ? (comparison[0] += 0.5 * temp1) : (comparison[0] += 0.5 * temp2);
				}
				for (int i = 0; i < NUM_USER; i++)
				{
					double temp1 = pow(H[i][j][0][sch] - estimate[(i + 1) % NUM_USER][j][sch][0], 2) + pow(H[i][j][1][sch] - estimate[(i + 1) % NUM_USER][j][sch][1], 2);
					double temp2 = pow(H[i][j][0][sch] + estimate[(i + 1) % NUM_USER][j][sch][0], 2) + pow(H[i][j][1][sch] + estimate[(i + 1) % NUM_USER][j][sch][1], 2);
					temp1 < temp2 ? (comparison[1] += 0.5 * temp1) : (comparison[1] += 0.5 * temp2);
				}
				if (comparison[0] > comparison[1])
				{
					Swaping(&estimate[0][j][sch], &estimate[1][j][sch]);
				}
			}
			else if (NUM_USER == 3)
			{
				double comparison[6] = { 0 };
				for (int k = 0; k < 6; k++)
				{
					for (int i = 0; i < NUM_USER; i++)
					{
						double temp1 = pow(H[i][j][0][sch] - estimate[i][j][sch][0], 2) + pow(H[i][j][1][sch] - estimate[i][j][sch][1], 2);
						double temp2 = pow(H[i][j][0][sch] + estimate[i][j][sch][0], 2) + pow(H[i][j][1][sch] + estimate[i][j][sch][1], 2);
						temp1 < temp2 ? (comparison[k] += 0.5 * temp1) : (comparison[k] += 0.5 * temp2);
					}
					Swaping(&estimate[k % NUM_USER][j][sch], &estimate[(k + 1) % NUM_USER][j][sch]);
				}
				int min = 0;
				for (int i = 1; i < 6; i++)
				{
					if (comparison[i] < comparison[min]) min = i;
				}
				for (int i = 0; i < min; i++)
				{
					Swaping(&estimate[i % NUM_USER][j][sch], &estimate[(i + 1) % NUM_USER][j][sch]);
				}
			}
		}
	}
	if (SLIDING) // adjustment for the same phase offset
	{
		for (int sch = 0; sch < FFT_POINT; sch++)
		{
			for (int i = 0; i < NUM_USER; i++)
			{
				for (int j = 1 + (WINDOW_SIZE - 1) / 2; j <= (FFT_SEGMENT + DIFF_ENC - WINDOW_SIZE) + (WINDOW_SIZE - 1) / 2; j++)
				{
					double temp1 = pow(estimate[i][j - 1][sch][0] - estimate[i][j][sch][0], 2) + pow(estimate[i][j - 1][sch][1] - estimate[i][j][sch][1], 2);
					double temp2 = pow(estimate[i][j - 1][sch][0] + estimate[i][j][sch][0], 2) + pow(estimate[i][j - 1][sch][1] + estimate[i][j][sch][1], 2);
					if (temp1 > temp2)
					{
						estimate[i][j][sch][0] *= -1;
						estimate[i][j][sch][1] *= -1;
					}
				}
				for (int j = 0; j < (WINDOW_SIZE - 1) / 2; j++)
				{
					estimate[i][j][sch][0] = estimate[i][(WINDOW_SIZE - 1) / 2][sch][0];
					estimate[i][j][sch][1] = estimate[i][(WINDOW_SIZE - 1) / 2][sch][1];
				}
				for (int j = FFT_SEGMENT + DIFF_ENC - WINDOW_SIZE + (WINDOW_SIZE - 1) / 2 + 1; j < FFT_SEGMENT + DIFF_ENC; j++)
				{
					estimate[i][j][sch][0] = estimate[i][FFT_SEGMENT + DIFF_ENC - WINDOW_SIZE + (WINDOW_SIZE - 1) / 2][sch][0];
					estimate[i][j][sch][1] = estimate[i][FFT_SEGMENT + DIFF_ENC - WINDOW_SIZE + (WINDOW_SIZE - 1) / 2][sch][1];
				}
			}
		}
	}
}

bool ConditionCheck(double **centroid, int *groupSize, int &count)
{
	for (int i = 0; i < GROUP_SIZE; i++)
	{
		if (groupSize[i] == 0) return true; // empty cluster
	}
	double average = 0;
	for (int i = 0; i < 2; i++)
	{
		double temp = 0;
		for (int j = 0; j < GROUP_SIZE; j++)
		{
			temp += centroid[j][i];
		}
		average += pow(temp / GROUP_SIZE, 2);
	}
	if (sqrt(average) > 0.05 + count*0.001) // ad-hoc
	{
		count++;
		return true;
	}
	return false;
}

void EMClustering(double ***postRx, double **centroid, double **softAssign, double variance, int sch, int offset)
{
	int it = 0;
	while (it++ < 1000)
	{
		if (SLIDING)
		{
			int i = 0;
			for (int k = 0; k < 2 * WINDOW_SIZE; k++)
			{
				if (k < WINDOW_SIZE)
					i = k + offset;
				else
					i = k + offset + FFT_SEGMENT + DIFF_ENC - WINDOW_SIZE;
		
				double sum = 0;
				for (int k = 0; k < GROUP_SIZE; k++)
				{
					softAssign[i][k] = exp(-pow(postRx[i][0][sch] - centroid[k][0], 2) / (2. * variance)) * exp(-pow(postRx[i][1][sch] - centroid[k][1], 2) / (2. * variance));
					if (softAssign[i][k] < 1e-300) softAssign[i][k] = 1e-300;
					sum += softAssign[i][k];
				}
				for (int k = 0; k < GROUP_SIZE; k++)
				{
					softAssign[i][k] /= sum;
					if (softAssign[i][k] < 1e-300) softAssign[i][k] = 1e-300;
				}
			}
		}
		else
		{
			for (int i = offset; i < 2*(FFT_SEGMENT + DIFF_ENC); i++)
			{
				double sum = 0;
				for (int k = 0; k < GROUP_SIZE; k++)
				{
					softAssign[i][k] = exp(-pow(postRx[i][0][sch] - centroid[k][0], 2) / (2. * variance)) * exp(-pow(postRx[i][1][sch] - centroid[k][1], 2) / (2. * variance));
					if (softAssign[i][k] < 1e-300) softAssign[i][k] = 1e-300;
					sum += softAssign[i][k];
				}
				for (int k = 0; k < GROUP_SIZE; k++)
				{
					softAssign[i][k] /= sum;
					if (softAssign[i][k] < 1e-300) softAssign[i][k] = 1e-300;
				}
			}
		}
		
		double variation = 0;
		for (int k = 0; k < GROUP_SIZE; k++)
		{
			double eMean[2] = { 0 }, nFactor = 0;
			for (int i = offset*SLIDING; i < (SLIDING ? offset + WINDOW_SIZE : FFT_SEGMENT + DIFF_ENC); i++)
			{
				eMean[0] += softAssign[i][k] * postRx[i][0][sch];
				eMean[1] += softAssign[i][k] * postRx[i][1][sch];
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

void printchannel(double ****h)
{
	FILE* ch_x = fopen("ch_x.txt", "w");
	FILE* ch_y = fopen("ch_y.txt", "w");

	for (int i = 0; i < FFT_SEGMENT + DIFF_ENC;i++)
	{
		fprintf(ch_x, "%e ", h[0][i][0][0]);
		fprintf(ch_y, "%e ", h[0][i][0][1]);
	}

	fclose(ch_x);
	fclose(ch_y);
	system("pause");
}