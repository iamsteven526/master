#define _USE_MATH_DEFINES
#include <iostream>
#include <random>
#include <fstream>
#include <algorithm>
#include <unordered_set>
#include "parameters.h"
#include <cstring>

#pragma warning(disable:4996)
using namespace std;

void KmeansClustering(double** rx, double** centroid, int** group, int* groupSize, double* distList, double* variation, double** softAssign, double variance, double** estimate, long double& itCount, double** chCoef)
{
	bool reset = true;
	int count = 0, it = 0;
//	double** sample =  new double*[2 * Cluster_Len];
//	for (int i = 0; i < 2* Cluster_Len; i++)
//		sample[i] = new double[2];

	double** sample =  new double*[Cluster_Len];
	for (int i = 0; i < Cluster_Len; i++)
		sample[i] = new double[2];

	for (int i = 0; i < Cluster_Len/2; i++)
	{
		sample[i][0] = rx[i][0];
		sample[i][1] = rx[i][1];
		sample[i+ Cluster_Len / 2][0] = -rx[i][0];
		sample[i+ Cluster_Len / 2][1] = -rx[i][1];
//		sample[i + Cluster_Len][0] = -rx[i][0];
//		sample[i + Cluster_Len][1] = -rx[i][1];
	}

	
	//system("pause");

	while (true)
	{
		itCount++;
		//if (count > 500)
		//	break;
		if (it++ > 10000)
		{
			for (int i = 0; i < NUM_USER; i++)
			{
				estimate[i][0] = estimate[i][1] = 0;
			}
			return;
		}
		if (reset)
		{
			InitialSeeding(sample, centroid, group, groupSize, distList);
			reset = false;
		}
		Grouping(sample, centroid, group, groupSize);
		CentroidRenewal(sample, centroid, group, groupSize, variation);
		bool convergence = true;
		for (int i = 0; i < GROUP_SIZE; i++)
		{
			if (variation[i] > 1e-5)
			{
				convergence = false;
				break;
			}
		}
		/*if (convergence)
		{
			//reset = ConditionCheck(centroid, groupSize, count);
			reset = ConditionCheck1(centroid, groupSize, variance, rx, estimate, chCoef, count, group);
			if (!reset)
			{
				
				if (EM_GMM) EMClustering(sample, centroid, softAssign, variance,chCoef,estimate);
				
				CoefEstimation(centroid, estimate, sample, chCoef);
				return;
			}

			//	printdata(rx, centroid,chCoef);
		}*/
		if (convergence == 1 && Check_initial_modified(sample, centroid, group, groupSize, variance, estimate, chCoef) == 1)
		{
			//Check_initial(rx, centroid, group, groupSize,variance);
			//reset = ConditionCheck1(centroid, groupSize,variance, rx, estimate, chCoef);
			reset = ConditionCheck(centroid, groupSize, count);
			if (!reset)
			{
				//if (EM_GMM) EMClustering(rx, centroid, softAssign, variance);
				CoefEstimation(centroid, estimate,sample,chCoef);
				return;
			}

		}
	}
	for (int i = 0; i < Cluster_Len; i++)
		delete[]sample[i];
	delete[] sample;
}

void InitialSeeding(double** sample, double** centroid, int** group, int* groupSize, double* distList)
{
	random_device seed;
	mt19937 generator(seed());
	uniform_real_distribution<double> uniform(0, 1);
	if (INI_METHOD == 1 || INI_METHOD == 2) // k-means++
	{
		uniform_int_distribution<int> ranPick(0, BLOCK_LEN - 1);
		int iniPick = ranPick(generator);
		centroid[0][0] = sample[iniPick][0]; centroid[0][1] = sample[iniPick][1];
		for (int i = 1; i < GROUP_SIZE; i++)
		{
			double sum = 0;
			for (int j = 0; j < BLOCK_LEN; j++) // distance from the closest centroid
			{
				distList[j] = EuclideanDistance(sample[j], centroid[0]);
				for (int k = 1; k < i; k++)
				{
					double dist = EuclideanDistance(sample[j], centroid[k]);
					if (dist < distList[j]) distList[j] = dist;
				}
				distList[j] = pow(distList[j], 2);
				sum += distList[j];
			}
			if (INI_METHOD == 2)
			{
				double list[2][BLOCK_LEN];
				memcpy(list[1], distList, sizeof(double) * BLOCK_LEN);
				for (int j = 0; j < BLOCK_LEN; j++)
				{
					list[0][j] = j;
				}
				for (int j = 0; j < BLOCK_LEN - 1; j++)// sort
				{
					for (int k = BLOCK_LEN - 1; k > j; k--)
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
				for (int j = 0; j < BLOCK_LEN - (BLOCK_LEN / GROUP_SIZE); j++)
				{
					sum -= distList[(int)list[0][j]];
					distList[(int)list[0][j]] = 0;
				}
			}
			double ranNum = uniform(generator) * sum;
			for (int j = 0; j < BLOCK_LEN; j++) // select with weighted probability
			{
				ranNum -= distList[j];
				if (ranNum <= 0)
				{
					centroid[i][0] = sample[j][0]; centroid[i][1] = sample[j][1];
					break;
				}
			}
		}
		//Check_initial(sample, centroid, group, groupSize);
	}
	else
	{
		printf("\nPARAMETER SETTING IS WRONG\n");
		system("pause");
	}
}

void Check_initial(double** rx, double** centroid, int** group, int* groupSize)
{

	bool tr = 0;
	for (int i = 0; i < GROUP_SIZE - 1; i++)
	{
		for (int j = i + 1; j < GROUP_SIZE; j++)
		{
			if (EuclideanDistance(centroid[i], centroid[j]) < 0.05)
			{
				tr = 1;
				break;
			}
		}
		if (tr == 1)
			break;
	}

	if (tr == 0)
		return;
	else
	{
	//	cout << 123;
	//	printinitial(rx, centroid);
		Grouping(rx, centroid, group, groupSize);

	/*		cout << endl << endl;
			for (int i=0;i<GROUP_SIZE;i++)
				cout << groupSize[i]<<" ";

			cout << endl;*/

		int min = 0, max = 0;
		int max_val = 0, min_val = 1e10;
		for (int i = 0; i < GROUP_SIZE; i++)
		{
			if (groupSize[i] < min_val)
			{
				min = i; min_val = groupSize[i];
			}
			if (groupSize[i] > max_val)
			{
				max = i; max_val = groupSize[i];
			}
		}

		//		cout << min_val << " " << max_val;

		if (groupSize[max] > groupSize[min] * 3)
		{
			centroid[min][0] = centroid[max][0];
			centroid[min][1] = centroid[max][1];
		}
		//		printinitial(rx, centroid);

	}
	//	system("pause");
}

void Grouping(double** sample, double** centroid, int** group, int* groupSize)
{
	memset(groupSize, 0, sizeof(int) * GROUP_SIZE);
	double dist[GROUP_SIZE];
	for (int i = 0; i < Cluster_Len; i++)
	{
		for (int j = 0; j < GROUP_SIZE; j++)
		{
			dist[j] = EuclideanDistance(sample[i], centroid[j]);
		}
		int min = 0;
		for (int j = 1; j < GROUP_SIZE; j++)
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

void CentroidRenewal(double** sample, double** centroid, int** group, int* groupSize, double* variation)
{
	memset(variation, 0, sizeof(double) * GROUP_SIZE);
	for (int i = 0; i < GROUP_SIZE; i++)
	{
		for (int j = 0; j < 2; j++) // real and imaginary
		{
			double temp = 0;
			for (int k = 0; k < groupSize[i]; k++)
			{
				temp += sample[group[i][k]][j];
			}
			temp /= groupSize[i];
			variation[i] += pow(centroid[i][j] - temp, 2);
			centroid[i][j] = temp;
		}
		variation[i] = sqrt(variation[i]);
	}
}

void CoefEstimation(double** centroid, double** estimate, double** sample, double** chCoef)
{
	if (NUM_USER == 1)
	{
		vector<vector<double>> testimate(2,vector<double> (2));
		
		double** phasor = new double* [GROUP_SIZE];
		for (int i = 0; i < GROUP_SIZE; i++)
		{
			phasor[i] = new double[2]; // amplitude and phase
		}
		for (int i = 0; i < GROUP_SIZE; i++) // estimation of superimposed signal
		{
			phasor[i][0] = sqrt(pow(centroid[i][0], 2) + pow(centroid[i][1], 2)); // amplitude
			phasor[i][1] = (centroid[i][1] > 0) ? (acos(centroid[i][0] / phasor[i][0])) : (2. * M_PI - acos(centroid[i][0] / phasor[i][0])); // phase
		}


		for (int i = 0; i < GROUP_SIZE - 1; i++) // clockwise sorting
		{
			for (int j = GROUP_SIZE - 1; j > i; j--)
			{
				if (phasor[j][1] < phasor[j - 1][1])
				{
					Swaping(&phasor[j], &phasor[j - 1]);
					Swaping(&centroid[j], &centroid[j - 1]);
				}
			}
		}

		double tphasor[2];// amplitude and phase
		tphasor[0] = 0.25*(phasor[0][0]+ phasor[1][0]+ phasor[2][0]+ phasor[3][0]); // amplitude
		tphasor[1] = 0.5*(phasor[0][1]+phasor[1][1]); // phase
		

		estimate[0][0] = tphasor[0] * cos(tphasor[1]);
		estimate[0][1] = tphasor[0] * sin(tphasor[1]);
	}
	else if (NUM_USER == 2)
	{
		int flag[GROUP_SIZE] = { 0 }, count = 0;
		int** pair = new int* [GROUP_SIZE >> 1];
		for (int i = 0; i < GROUP_SIZE >> 1; i++)
			pair[i] = new int[2];
		double** sum = new double* [GROUP_SIZE - 1];
		for (int i = 0; i < GROUP_SIZE - 1; i++)
		{
			sum[i] = new double[GROUP_SIZE - i - 1];
		}
		vector<vector<double>> centroid_permutation;


		for (int i = 0; i < GROUP_SIZE - 1; i++) // compute the sum of each pair //sum �`�@��7+6+5+4+3+2+1
		{
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
				for (int j = 0; j < GROUP_SIZE - i - 1; j++)
				{
					if (sum[i][j] < min && flag[i] == 0 && flag[i + j + 1] == 0) //�Ysum�p��min�h���@��pair
					{
						min = sum[i][j];
						pair[count][0] = i;
						pair[count][1] = i + j + 1;
					}
				}
			}
			flag[pair[count][0]] = flag[pair[count][1]] = 1;
			count++;
		}

		
		QPSKcheck(centroid, pair, centroid_permutation,estimate);
		for (int i = 0; i < GROUP_SIZE >> 1; i++)
			delete[] pair[i];

		for (int i = 0; i < GROUP_SIZE - 1; i++)
			delete[] sum[i];
		
		delete[] pair;
		delete[] sum;

	//	printdata(sample, estimate, chCoef, centroid);		
	}
	else if (NUM_USER == 3)
	{

	}
}

void MSEComparison(double** chCoef, double** estimate, double& mse, double** sample, double** centroid)
{

	if (NUM_USER == 1)
	{
		vector<double> tphasor(2);// amplitude and phase
		tphasor[0] = sqrt(pow(estimate[0][0], 2) + pow(estimate[0][1], 2)) ; // amplitude
		tphasor[1] = (estimate[0][1] > 0) ? (acos(estimate[0][0] / tphasor[0])) : (2. * M_PI - acos(estimate[0][0] / tphasor[0])); // phase
		

		vector<double> temp(4);
		temp[0] = pow(chCoef[0][0] - estimate[0][0], 2) + pow(chCoef[0][1] - estimate[0][1], 2);
		temp[1] = pow(chCoef[0][0] + estimate[0][0], 2) + pow(chCoef[0][1] + estimate[0][1], 2);
		temp[2] = pow(chCoef[0][0] - tphasor[0] * cos(tphasor[1] + M_PI / 2 ), 2) + pow(chCoef[0][1] - tphasor[0] * sin(tphasor[1] + M_PI / 2), 2);
		temp[3] = pow(chCoef[0][0] - tphasor[0] * cos(tphasor[1] + 3*M_PI / 2), 2) + pow(chCoef[0][1] - tphasor[0] * sin(tphasor[1] + 3 * M_PI / 2), 2);
		
		double reg,min=10;
		for (int i = 0; i < 4; i++)
		{
			if (temp[i] < min)
			{
				reg = i;
				min = temp[i];
			}
		}


		mse += 0.5*temp[reg];

		if (temp[reg] > 0.0005)
		{
			printdata(sample, estimate, chCoef, centroid);
		}
		tphasor.clear();
		temp.clear();
	}
	else if (NUM_USER == 2)
	{
		double temp;
		double reg=0;

		double min[2] = { 1e10,1e10 };
		//	Compare chCoef_user1 chCoef_user2 with estimate_user1
		for (int i = 0; i < NUM_USER; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				temp = pow(chCoef[i][2] * cos((M_PI / 2)*j + chCoef[i][4]) - estimate[0][0], 2) + pow(chCoef[i][2] * sin((M_PI / 2)*j + chCoef[i][4]) - estimate[0][1], 2);
				if (temp < min[i])
					min[i] = temp;
			}
		}

		// chCoef_user1 = estimate_user2, chCoef_user2 = estimate_user1
		if (min[1] < min[0])
		{
			reg += min[1];
			swap(estimate[0], estimate[1]);
			min[1] = 1e10;
			
			for (int j = 0; j < 4; j++)
			{
				temp = pow(chCoef[0][2] * cos((M_PI / 2)*j + chCoef[0][4]) - estimate[0][0], 2) + pow(chCoef[0][2] * sin((M_PI / 2)*j + chCoef[0][4]) - estimate[0][1], 2);
				if (temp < min[1])
					min[1] = temp;
			}
			reg += min[1];
		}
		// chCoef_user1 = estimate_user1, chCoef_user2 = estimate_user2
		else
		{
			reg += min[0];
			min[0] = 1e10;
			for (int j = 0; j < 4; j++)
			{
				temp = pow(chCoef[1][2] * cos((M_PI / 2) * j + chCoef[1][4]) - estimate[1][0], 2) + pow(chCoef[1][2] * sin((M_PI / 2) * j + chCoef[1][4]) - estimate[1][1], 2);
				if (temp < min[0])
					min[0] = temp;
			}
			reg += min[0];
		}

		mse += reg;

		//if(reg>0.01)
		//	printdata(sample, estimate, chCoef, centroid);	
	}
	else if (NUM_USER == 3)
	{
		
	}
}

bool ConditionCheck(double** centroid, int* groupSize, int& count)//true �~��, false ���U��
{
	for (int i = 0; i < GROUP_SIZE; i++)//�ˬdGROUP_SIZE��centroid�����I
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
	if (sqrt(average) > 0.05 + count * 0.001) // ad-hoc
	{
		count++;
		return true;
	}
	return false;
}

void EMClustering(double** sample, double** centroid, double** softAssign, double variance,double** chCoef,double** estimate)
{
	int it = 1;
	//cout << "=============================";
	while (it++ < 1000)
	{
		for (int i = 0; i < Cluster_Len; i++)
		{
			double sum = 0;
			for (int k = 0; k < GROUP_SIZE; k++)
			{
				softAssign[i][k] = exp(-pow(sample[i][0] - centroid[k][0], 2) / (2. * variance)) * exp(-pow(sample[i][1] - centroid[k][1], 2) / (2. * variance));
				if (softAssign[i][k] < NUMERIC_LIMIT) softAssign[i][k] = NUMERIC_LIMIT;
				sum += softAssign[i][k];
			}
			for (int k = 0; k < GROUP_SIZE; k++)
			{
				softAssign[i][k] /= sum;
				if (softAssign[i][k] < NUMERIC_LIMIT) softAssign[i][k] = NUMERIC_LIMIT;
			}
		}
		double diff = 0;
		for (int k = 0; k < GROUP_SIZE; k++)
		{
			double eMean[2] = { 0 }, nFactor = 0;
			for (int i = 0; i < Cluster_Len; i++)
			{
				eMean[0] += softAssign[i][k] * sample[i][0];
				eMean[1] += softAssign[i][k] * sample[i][1];
				nFactor += softAssign[i][k];
			}
			eMean[0] /= nFactor;
			eMean[1] /= nFactor;
			diff += pow(centroid[k][0] - eMean[0], 2) + pow(centroid[k][1] - eMean[1], 2);
			centroid[k][0] = eMean[0]; centroid[k][1] = eMean[1];

		//	cout << centroid[k][0] << " " << centroid[k][1]<<endl;

		//	if (centroid[k][0] < 0.05 && centroid[k][1] < 0.05)
		//		printdata(rx,estimate,chCoef,centroid);
		}
		//system("pause");
		if (diff < 1e-10) return;
	}
	
}

void printdata(double** sample, double** es, double** chCoef, double** centroid)
{
	FILE* rx_x = fopen("rx_x.txt", "w");
	FILE* rx_y = fopen("rx_y.txt", "w");
	FILE* ce_x = fopen("ce_x.txt", "w");
	FILE* ce_y = fopen("ce_y.txt", "w");
	FILE* ch_x = fopen("ch_x.txt", "w");
	FILE* ch_y = fopen("ch_y.txt", "w");

	FILE* centroid_x = fopen("centroid_x.txt", "w");
	FILE* centroid_y = fopen("centroid_y.txt", "w");
	

	for (int j = 0; j < Cluster_Len; j++)
	{
		fprintf(rx_x, "%e ", sample[j][0]);
	}

	for (int j = 0; j < Cluster_Len; j++)
	{
		fprintf(rx_y, "%e ", sample[j][1]);
	}


	for (int j = 0; j < GROUP_SIZE; j++)
	{
		fprintf(centroid_x, "%e ", centroid[j][0]);
	}

	for (int j = 0; j < GROUP_SIZE; j++)
	{
		fprintf(centroid_y, "%e ", centroid[j][1]);
	}

	vector<vector<double>> tphasor(NUM_USER,vector<double>(2));

	for (int i = 0; i < NUM_USER; i++)
	{
		tphasor[i][0] = sqrt(pow(es[i][0], 2) + pow(es[i][1], 2)); // amplitude
		tphasor[i][1] = (es[i][1] > 0) ? (acos(es[i][0] / tphasor[i][0])) : (2. * M_PI - acos(es[i][0] / tphasor[i][0])); // phase
	}
//	tphasor[1][0] = sqrt(pow(es[1][0], 2) + pow(es[1][1], 2)); // amplitude
//	tphasor[1][1] = (es[1][1] > 0) ? (acos(es[1][0] / tphasor[1][0])) : (2. * M_PI - acos(es[1][0] / tphasor[1][0])); // phase
	

	for (int j = 0; j < NUM_USER; j++)
	{
		fprintf(ce_x, "%e ", tphasor[j][0] * cos(tphasor[j][1] ));
		fprintf(ce_x, "%e ", tphasor[j][0] * cos(tphasor[j][1] + M_PI / 2));
		fprintf(ce_x, "%e ", tphasor[j][0] * cos(tphasor[j][1] + M_PI));
		fprintf(ce_x, "%e ", tphasor[j][0] * cos(tphasor[j][1] + 3*M_PI / 2));
	}



	for (int j = 0; j < NUM_USER; j++)
	{
		fprintf(ce_y, "%e ", tphasor[j][0] * sin(tphasor[j][1]));
		fprintf(ce_y, "%e ", tphasor[j][0] * sin(tphasor[j][1] + M_PI / 2));
		fprintf(ce_y, "%e ", tphasor[j][0] * sin(tphasor[j][1] + M_PI));
		fprintf(ce_y, "%e ", tphasor[j][0] * sin(tphasor[j][1] + 3*M_PI / 2));
	}


	for (int j = 0; j < NUM_USER; j++)
	{
		fprintf(ch_x, "%e ", chCoef[j][0]);
	}

	for (int j = 0; j < NUM_USER; j++)
	{
		fprintf(ch_y, "%e ", chCoef[j][1]);
	}

	fclose(rx_x);
	fclose(rx_y);
	fclose(ce_x);
	fclose(ce_y);
	fclose(ch_x);
	fclose(ch_y);
	fclose(centroid_x);
	fclose(centroid_y);
	

	system("pause");
}

void printinitial(double** sample, double** centroid)
{
	FILE* rx_x_intitial = fopen("rx_x_intitial.txt", "w");
	FILE* rx_y_intitial = fopen("rx_y_intitial.txt", "w");
	FILE* centroid_x_intitial = fopen("centroid_x_intitial.txt", "w");
	FILE* centroid_y_intitial = fopen("centroid_y_intitial.txt", "w");

	for (int j = 0; j < Cluster_Len; j++)
	{
		fprintf(rx_x_intitial, "%e ", sample[j][0]);
	}

	for (int j = 0; j < Cluster_Len; j++)
	{
		fprintf(rx_y_intitial, "%e ", sample[j] [1] );
	}


	for (int j = 0; j < GROUP_SIZE; j++)
	{
		if (centroid[j][0] == 0)
			continue;
		fprintf(centroid_x_intitial, "%e ", centroid[j][0]);
	}

	for (int j = 0; j < GROUP_SIZE; j++)
	{
		if (centroid[j][1] == 0)
			continue;
		fprintf(centroid_y_intitial, "%e ", centroid[j][1]);
	}

	fclose(rx_x_intitial);
	fclose(rx_y_intitial);
	fclose(centroid_x_intitial);
	fclose(centroid_y_intitial);

	system("pause");
}

void CentroidMSEComparison(double** chCoef, double** estimate, double& centroid_mse, double** sample, double** centroid)
{
	if (NUM_USER == 1)
	{
		vector<vector<double>> level_val(NUM_LEVEL, vector<double>(2));
		vector<double> tphasor(2);// amplitude and phase
		tphasor[0] = sqrt(pow(chCoef[0][0], 2) + pow(chCoef[0][1], 2)); // amplitude
		tphasor[1] = (chCoef[0][1] > 0) ? (acos(chCoef[0][0] / tphasor[0])) : (2. * M_PI - acos(chCoef[0][0] / tphasor[0])); // phase

		level_val[0][0] = tphasor[0] * cos(tphasor[1] + M_PI / 4);
		level_val[0][1] = tphasor[0] * sin(tphasor[1] + M_PI / 4);
		level_val[1][0] = -tphasor[0] * cos(tphasor[1] + M_PI / 4);
		level_val[1][1] = -tphasor[0] * sin(tphasor[1] + M_PI / 4 );
		level_val[2][0] = tphasor[0] * cos(tphasor[1] + 3 * M_PI / 4);
		level_val[2][1] = tphasor[0] * sin(tphasor[1] + 3 * M_PI / 4);
		level_val[3][0] = -tphasor[0] * cos(tphasor[1] + 3 * M_PI / 4);
		level_val[3][1] = -tphasor[0] * sin(tphasor[1] + 3 * M_PI / 4);

	
		double temp = 100, temp2;

		int min = 0;

		for (int i = 0; i < 4; i++)
		{
			temp = 100;
			min = i;
			for (int j = i; j < 4; j++)
			{
				double temp1 = pow(level_val[j][0] - centroid[i][0], 2) + pow(level_val[j][1] - centroid[i][1], 2);
				if (0.5 * temp1 < temp)
				{
					temp = 0.5 * temp1;
					min = j;
				}
			}


			centroid_mse += temp;

		//	if(temp>0.1)
		//		printdata(rx, estimate, chCoef, centroid);

			if (min != i)
			{
				temp2 = level_val[i][0];
				level_val[i][0] = level_val[min][0];
				level_val[min][0] = temp2;

				temp2 = level_val[i][1];
				level_val[i][1] = level_val[min][1];
				level_val[min][1] = temp2;
			}

		}
	}
	else if (NUM_USER == 2)
	{
		vector<vector<double>> level_val(NUM_LEVEL, vector<double>(2));


		level_val[0][0] = chCoef[1][2] * cos(chCoef[1][4] + M_PI / 4) + chCoef[0][2] * cos(chCoef[0][4] + M_PI / 4);
		level_val[0][1] = chCoef[1][2] * sin(chCoef[1][4] + M_PI / 4) + chCoef[0][2] * sin(chCoef[0][4] + M_PI / 4);
		level_val[1][0] = chCoef[1][2] * cos(chCoef[1][4] + M_PI / 4) - chCoef[0][2] * cos(chCoef[0][4] + M_PI / 4);
		level_val[1][1] = chCoef[1][2] * sin(chCoef[1][4] + M_PI / 4) - chCoef[0][2] * sin(chCoef[0][4] + M_PI / 4);
		level_val[2][0] = chCoef[1][2] * cos(chCoef[1][4] + M_PI / 4) + chCoef[0][2] * cos(chCoef[0][4] + 3 * M_PI / 4);
		level_val[2][1] = chCoef[1][2] * sin(chCoef[1][4] + M_PI / 4) + chCoef[0][2] * sin(chCoef[0][4] + 3 * M_PI / 4);
		level_val[3][0] = chCoef[1][2] * cos(chCoef[1][4] + M_PI / 4) - chCoef[0][2] * cos(chCoef[0][4] + 3 * M_PI / 4);
		level_val[3][1] = chCoef[1][2] * sin(chCoef[1][4] + M_PI / 4) - chCoef[0][2] * sin(chCoef[0][4] + 3 * M_PI / 4);

		level_val[4][0] = chCoef[1][2] * cos(chCoef[1][4] + 3* M_PI / 4) + chCoef[0][2] * cos(chCoef[0][4] + M_PI / 4);
		level_val[4][1] = chCoef[1][2] * sin(chCoef[1][4] + 3 * M_PI / 4) + chCoef[0][2] * sin(chCoef[0][4] + M_PI / 4);
		level_val[5][0] = chCoef[1][2] * cos(chCoef[1][4] + 3 * M_PI / 4) - chCoef[0][2] * cos(chCoef[0][4] + M_PI / 4);
		level_val[5][1] = chCoef[1][2] * sin(chCoef[1][4] + 3 * M_PI / 4) - chCoef[0][2] * sin(chCoef[0][4] + M_PI / 4);
		level_val[6][0] = chCoef[1][2] * cos(chCoef[1][4] + 3 * M_PI / 4) + chCoef[0][2] * cos(chCoef[0][4] + 3 * M_PI / 4);
		level_val[6][1] = chCoef[1][2] * sin(chCoef[1][4] + 3 * M_PI / 4) + chCoef[0][2] * sin(chCoef[0][4] + 3 * M_PI / 4);
		level_val[7][0] = chCoef[1][2] * cos(chCoef[1][4] + 3 * M_PI / 4) - chCoef[0][2] * cos(chCoef[0][4] + 3 * M_PI / 4);
		level_val[7][1] = chCoef[1][2] * sin(chCoef[1][4] + 3 * M_PI / 4) - chCoef[0][2] * sin(chCoef[0][4] + 3 * M_PI / 4);
		
		level_val[8][0] = -level_val[0][0];
		level_val[8][1] = -level_val[0][1];
		level_val[9][0] = -level_val[1][0];
		level_val[9][1] = -level_val[1][1];
		level_val[10][0] = -level_val[2][0];
		level_val[10][1] = -level_val[2][1];
		level_val[11][0] = -level_val[3][0];
		level_val[11][1] = -level_val[3][1];

		level_val[12][0] = -level_val[4][0];
		level_val[12][1] = -level_val[4][1];
		level_val[13][0] = -level_val[5][0];
		level_val[13][1] = -level_val[5][1];
		level_val[14][0] = -level_val[6][0];
		level_val[14][1] = -level_val[6][1];
		level_val[15][0] = -level_val[7][0];
		level_val[15][1] = -level_val[7][1];
		

		double temp = 100, temp2;

		int min = 0;

		for (int i = 0; i < 16; i++)
		{
			temp = 100;
			min = i;
			for (int j = i; j < 16; j++)
			{
				double temp1 = pow(level_val[j][0] - centroid[i][0], 2) + pow(level_val[j][1] - centroid[i][1], 2);
				if (0.5 * temp1 < temp)
				{
					temp = 0.5 * temp1;
					min = j;
				}
			}

			centroid_mse += temp;

			//if (temp > 0.1)
			//	printdata(sample, estimate, chCoef, centroid);

			if (min != i)
			{
				temp2 = level_val[i][0];
				level_val[i][0] = level_val[min][0];
				level_val[min][0] = temp2;

				temp2 = level_val[i][1];
				level_val[i][1] = level_val[min][1];
				level_val[min][1] = temp2;
			}

		}
	//	printdata(rx, estimate, chCoef, centroid);
	}
	else if (NUM_USER == 3)
	{

	}
}

void QPSKcheck(double** centroid,int **pair,vector<vector<double>> & centroid_per, double  **estimate)
{	
	// Generate first and second pair square sequence
	Generate_square_sequence(centroid_per,pair,centroid,8,2);

	// Record first Square-Sequence number to generate third Square-Sequence 
	unordered_set<int> s;
	int ** new_pair=new int*[4];
	for (int i = 0; i < 4; i++)
		new_pair[i] = new int[2];
	

	for (int j = 0; j < 4; j++)
	{
		s.insert(centroid_per[0][j]);
		s.insert(centroid_per[1][j]);
	}

	int c = 0;
	for (int i = 0; i < 8; i++)
	{
		if (s.count(pair[i][0]) == 0 && s.count(pair[i][1]) == 0)
		{
			new_pair[c][0] = pair[i][0];
			new_pair[c][1] = pair[i][1];
			c++;
		}
	}

	Generate_square_sequence(centroid_per, new_pair, centroid,4,1);
	
	// Record second Square-Sequence number to generate fourth Square-Sequence 
	s.clear();
	


	for (int j = 0; j < 4; j++)
	{
		s.insert(centroid_per[2][j]);
		s.insert(centroid_per[3][j]);
	}

	c = 0;
	for (int i = 0; i < 8; i++)
	{
		if (s.count(pair[i][0]) == 0 && s.count(pair[i][1]) == 0)
		{
			new_pair[c][0] = pair[i][0];
			new_pair[c][1] = pair[i][1];
			c++;
		}
	}

	Generate_square_sequence(centroid_per, new_pair, centroid, 4, 1);


	// Calculate estimation of centroid by Square-Sequence
	// The sequence order of Square-Sequence : the first and third sequence are same user, the second and fourth sequence are same user.
	
	vector<vector<double>> testimate(4,vector<double>(2));
	// User A-1 : Average the same phase( pi ) sequence
	testimate[0][0] = (centroid[int(centroid_per[0][0])][0] + centroid[int(centroid_per[0][1])][0] + centroid[int(centroid_per[0][2])][0] + centroid[int(centroid_per[0][3])][0]) / 4;  
	testimate[0][1] = (centroid[int(centroid_per[0][0])][1] + centroid[int(centroid_per[0][1])][1] + centroid[int(centroid_per[0][2])][1] + centroid[int(centroid_per[0][3])][1]) / 4;  
	testimate[0][0] = 0.5 * testimate[0][0] - (centroid[int(centroid_per[1][0])][0] + centroid[int(centroid_per[1][1])][0] + centroid[int(centroid_per[1][2])][0] + centroid[int(centroid_per[1][3])][0]) / 8;  
	testimate[0][1] = 0.5 * testimate[0][1] - (centroid[int(centroid_per[1][0])][1] + centroid[int(centroid_per[1][1])][1] + centroid[int(centroid_per[1][2])][1] + centroid[int(centroid_per[1][3])][1]) / 8;  
	
	// User B-1 : 
	testimate[1][0] = (centroid[int(centroid_per[2][0])][0] + centroid[int(centroid_per[2][1])][0] + centroid[int(centroid_per[2][2])][0] + centroid[int(centroid_per[2][3])][0]) / 4;  
	testimate[1][1] = (centroid[int(centroid_per[2][0])][1] + centroid[int(centroid_per[2][1])][1] + centroid[int(centroid_per[2][2])][1] + centroid[int(centroid_per[2][3])][1]) / 4;  
	testimate[1][0] = 0.5 * testimate[1][0] - (centroid[int(centroid_per[3][0])][0] + centroid[int(centroid_per[3][1])][0] + centroid[int(centroid_per[3][2])][0] + centroid[int(centroid_per[3][3])][0]) / 8;  
	testimate[1][1] = 0.5 * testimate[1][1] - (centroid[int(centroid_per[3][0])][1] + centroid[int(centroid_per[3][1])][1] + centroid[int(centroid_per[3][2])][1] + centroid[int(centroid_per[3][3])][1]) / 8;  

	// User A-2 : (pi/2) to User A-1 
	testimate[2][0] = (centroid[int(centroid_per[4][0])][0] + centroid[int(centroid_per[4][1])][0] + centroid[int(centroid_per[4][2])][0] + centroid[int(centroid_per[4][3])][0]) / 4;
	testimate[2][1] = (centroid[int(centroid_per[4][0])][1] + centroid[int(centroid_per[4][1])][1] + centroid[int(centroid_per[4][2])][1] + centroid[int(centroid_per[4][3])][1]) / 4;
	testimate[2][0] = 0.5 * testimate[2][0] - (centroid[int(centroid_per[5][0])][0] + centroid[int(centroid_per[5][1])][0] + centroid[int(centroid_per[5][2])][0] + centroid[int(centroid_per[5][3])][0]) / 8;
	testimate[2][1] = 0.5 * testimate[2][1] - (centroid[int(centroid_per[5][0])][1] + centroid[int(centroid_per[5][1])][1] + centroid[int(centroid_per[5][2])][1] + centroid[int(centroid_per[5][3])][1]) / 8;

	// User B-2 : (pi/2) to User B-1
	testimate[3][0] = (centroid[int(centroid_per[6][0])][0] + centroid[int(centroid_per[6][1])][0] + centroid[int(centroid_per[6][2])][0] + centroid[int(centroid_per[6][3])][0]) / 4;
	testimate[3][1] = (centroid[int(centroid_per[6][0])][1] + centroid[int(centroid_per[6][1])][1] + centroid[int(centroid_per[6][2])][1] + centroid[int(centroid_per[6][3])][1]) / 4;
	testimate[3][0] = 0.5 * testimate[3][0] - (centroid[int(centroid_per[7][0])][0] + centroid[int(centroid_per[7][1])][0] + centroid[int(centroid_per[7][2])][0] + centroid[int(centroid_per[7][3])][0]) / 8;
	testimate[3][1] = 0.5 * testimate[3][1] - (centroid[int(centroid_per[7][0])][1] + centroid[int(centroid_per[7][1])][1] + centroid[int(centroid_per[7][2])][1] + centroid[int(centroid_per[7][3])][1]) / 8;


	vector<vector<double>> phasor(4,vector<double> (2));
	vector<vector<double>> tphasor(2, vector<double>(2));

	for (int i = 0; i < 4; i++) // estimation of superimposed signal
	{
		phasor[i][0] = sqrt(pow(testimate[i][0], 2) + pow(testimate[i][1], 2)); // amplitude
		phasor[i][1] = (testimate[i][1] > 0) ? (acos(testimate[i][0] / phasor[i][0])) : (2. * M_PI - acos(testimate[i][0] / phasor[i][0])); // phase
	}
	
	tphasor[0][0] = 0.5 * (phasor[0][0] + phasor[2][0]);
	tphasor[0][1] = 0.5 * (phasor[0][1] + phasor[2][1]);

	tphasor[1][0] = 0.5 * (phasor[1][0] + phasor[3][0]);
	tphasor[1][1] = 0.5 * (phasor[1][1] + phasor[3][1]);

	estimate[0][0] = tphasor[0][0] * cos(tphasor[0][1]);
	estimate[0][1] = tphasor[0][0] * sin(tphasor[0][1]);

	estimate[1][0] = tphasor[1][0] * cos(tphasor[1][1]);
	estimate[1][1] = tphasor[1][0] * sin(tphasor[1][1]);


	for (int i = 0; i < 4; i++)
	{
		delete[] new_pair[i];
	}

	delete[] new_pair;

	phasor.clear();
	tphasor.clear();
}

void Generate_square_sequence(vector<vector<double>> & centroid_per, int ** pair, double **centroid,int pair_size,int sequence_num)
{
	double min = 1e10;
	double bound = 0.01;

	vector<vector<int>> sorted_diff_1th; // the number of ���y�T����
	vector<vector<int>> sorted_diff_2th; // the number of �����
	vector<vector<double>> temp_centroid_per;

	while (sorted_diff_2th.size() < sequence_num)  // check�Ĥ@�ӵ��y�T����
	{
		bound += 0.01;
		sorted_diff_1th.clear();
		sorted_diff_2th.clear();

		for (int i = 1; i < pair_size-1; i++)
		{
			for (int j = i + 1; j < pair_size; j++)
			{
				vector<double> sum_compare(4); // four state of 00 01 10 11
				sum_compare[0] = abs(EuclideanDistance(centroid[pair[0][0]], centroid[pair[i][0]]) - EuclideanDistance(centroid[pair[0][0]], centroid[pair[j][0]]));
				sum_compare[1] = abs(EuclideanDistance(centroid[pair[0][0]], centroid[pair[i][0]]) - EuclideanDistance(centroid[pair[0][0]], centroid[pair[j][1]]));
				sum_compare[2] = abs(EuclideanDistance(centroid[pair[0][0]], centroid[pair[i][1]]) - EuclideanDistance(centroid[pair[0][0]], centroid[pair[j][0]]));
				sum_compare[3] = abs(EuclideanDistance(centroid[pair[0][0]], centroid[pair[i][1]]) - EuclideanDistance(centroid[pair[0][0]], centroid[pair[j][1]]));

				for (int k = 0; k < 4; k++)
				{
					if (sum_compare[k] < bound)
					{
						vector<int> temp(3);
						temp[0] = i; temp[1] = j; temp[2] = k;
						sorted_diff_1th.push_back(temp);
						temp.clear();
					}
				}
				sum_compare.clear();
			}
		}

	/*	for (int i = 0; i < sorted_diff_1th.size(); i++)
		{
			for (int j = 0; j < 3; j++)
			{
				cout << sorted_diff_1th[i][j] << " ";
			}
			cout << endl;
		}*/


		for (int i = 0; i < sorted_diff_1th.size(); i++) // check�ĤG�ӵ��y�T����
		{
			for (int j = 1; j < pair_size; j++)
			{
				if (j == sorted_diff_1th[i][0] || j == sorted_diff_1th[i][1])
					continue;
				else
					for (int k = 0; k < 2; k++)
					{
						double reg;
						if (sorted_diff_1th[i][2] == 0)
							reg = abs(EuclideanDistance(centroid[pair[sorted_diff_1th[i][0]][0]], centroid[pair[j][k]]) - EuclideanDistance(centroid[pair[sorted_diff_1th[i][1]][0]], centroid[pair[j][k]]));

						else if (sorted_diff_1th[i][2] == 1)
							reg = abs(EuclideanDistance(centroid[pair[sorted_diff_1th[i][0]][0]], centroid[pair[j][k]]) - EuclideanDistance(centroid[pair[sorted_diff_1th[i][1]][1]], centroid[pair[j][k]]));

						else if (sorted_diff_1th[i][2] == 2)
							reg = abs(EuclideanDistance(centroid[pair[sorted_diff_1th[i][0]][1]], centroid[pair[j][k]]) - EuclideanDistance(centroid[pair[sorted_diff_1th[i][1]][0]], centroid[pair[j][k]]));

						else
							reg = abs(EuclideanDistance(centroid[pair[sorted_diff_1th[i][0]][1]], centroid[pair[j][k]]) - EuclideanDistance(centroid[pair[sorted_diff_1th[i][1]][1]], centroid[pair[j][k]]));


						if (reg < bound)
						{
							vector<int> temp(5);
							temp[0] = sorted_diff_1th[i][0]; temp[1] = sorted_diff_1th[i][1]; temp[2] = sorted_diff_1th[i][2];
							temp[3] = j; temp[4] = k;
							sorted_diff_2th.push_back(temp);
							temp.clear();
						}
					}
			}
		}
	}

	/*for (int i = 0; i < sorted_diff_2th.size(); i++)
	{
		for (int j = 0; j < 5; j++)
		{
			cout << sorted_diff_2th[i][j] << " ";
		}
		cout << endl;
	}*/

	bound = 0.001;
	int count = 0;
	while (temp_centroid_per.size() < sequence_num)//�ഫ����νs�� and check�����
	{
		temp_centroid_per.clear();
		for (int i = 0; i < sorted_diff_2th.size(); i++)
		{
			vector<int> temp(4);
			temp[0] = pair[0][0];

			if (sorted_diff_2th[i][2] == 0)
			{
				temp[1] = pair[sorted_diff_2th[i][0]][0];
				temp[3] = pair[sorted_diff_2th[i][1]][0];
			}
			else if (sorted_diff_2th[i][2] == 1)
			{
				temp[1] = pair[sorted_diff_2th[i][0]][0];
				temp[3] = pair[sorted_diff_2th[i][1]][1];
			}
			else if (sorted_diff_2th[i][2] == 2)
			{
				temp[1] = pair[sorted_diff_2th[i][0]][1];
				temp[3] = pair[sorted_diff_2th[i][1]][0];
			}
			else
			{
				temp[1] = pair[sorted_diff_2th[i][0]][1];
				temp[3] = pair[sorted_diff_2th[i][1]][1];
			}

			if (sorted_diff_2th[i][4] == 0)
				temp[2] = pair[sorted_diff_2th[i][3]][0];
			else
				temp[2] = pair[sorted_diff_2th[i][3]][1];


			double reg, sum = 0;

			reg = abs(EuclideanDistance(centroid[temp[0]], centroid[temp[1]]) - EuclideanDistance(centroid[temp[1]], centroid[temp[2]]));
			sum += reg;
			if (reg > bound)
				continue;

			reg = abs(EuclideanDistance(centroid[temp[1]], centroid[temp[2]]) - EuclideanDistance(centroid[temp[2]], centroid[temp[3]]));
			sum += reg;
			if (reg > bound)
				continue;

			reg = abs(EuclideanDistance(centroid[temp[2]], centroid[temp[3]]) - EuclideanDistance(centroid[temp[3]], centroid[temp[0]]));
			sum += reg;
			if (reg > bound)
				continue;

			reg = abs(EuclideanDistance(centroid[temp[0]], centroid[temp[2]]) - EuclideanDistance(centroid[temp[1]], centroid[temp[3]]));
			sum += reg;
			if (reg > bound)
				continue;

			vector<double> reg2; // square �W�۹諸�s��

			for (int i = 0; i < 4; i++)
				reg2.push_back(double(temp[i]));

			reg2.push_back(sum);

			temp_centroid_per.push_back(reg2);
		}
		bound += 0.001;

		if (temp_centroid_per.size() > sequence_num)
		{
			/*for (int i = 0; i < temp_centroid_per.size(); i++)
			{
				for (int j = 0; j < 5; j++)
				{
					cout << temp_centroid_per[i][j] << " ";
				}
				cout << endl;
			}*/
			//cout << "-----------" << endl;
			for (int i = 0; i < temp_centroid_per.size()-1; i++)
			{
				for (int j = i+1; j < temp_centroid_per.size(); j++)
				{
					if (temp_centroid_per[i][4] > temp_centroid_per[j][4])
					{
						for (int k = 0; k < 5; k++)
						{
							double reg;
							reg = temp_centroid_per[i][k];
							temp_centroid_per[i][k] = temp_centroid_per[j][k];
							temp_centroid_per[j][k] = reg;
						}
					}
				}
			}
			/*for (int i = 0; i < temp_centroid_per.size(); i++)
			{
				for (int j = 0; j < 5; j++)
				{
					cout << temp_centroid_per[i][j] << " ";
				}
				cout << endl;
			}*/
			//cout << "---------------------" << endl;
		}
	}

	vector<int> pair_table( NUM_LEVEL);
	for (int i = 0; i < NUM_LEVEL; i++)
	{
		for (int j = 0; j < pair_size; j++)
		{
			if (pair[j][0] == i)
				pair_table[i] = pair[j][1];
			else if (pair[j][1] == i)
				pair_table[i] = pair[j][0];
			else
				continue;
		}
	}

//	for (auto val : pair_table)
//		cout << val << " ";
	//system("pause");
	for (int i = 0; i < sequence_num; i++) // generate symmetry square
	{
		vector<double> reg;
		for (int j = 0; j < 4; j++)
			reg.push_back(temp_centroid_per[i][j]);
		centroid_per.push_back(reg);
		reg.clear();

		for (int j = 0; j < 4; j++)
			reg.push_back(pair_table[temp_centroid_per[i][j]]);

		centroid_per.push_back(reg);
	}
}

bool Check_initial_modified(double** rx, double** centroid, int** group, int* groupSize, double variance, double** estimate, double** chCoef)
{

	vector<double> mse(GROUP_SIZE);

	for (int i = 0; i < GROUP_SIZE; i++)
	{
		for (int j = 0; j < groupSize[i]; j++)
		{
			mse[i] += 0.5 * (pow(centroid[i][0] - rx[group[i][j]][0], 2) + pow(centroid[i][1] - rx[group[i][j]][1], 2));
		}
		mse[i] /= groupSize[i];
	}
	/*	double sum = 0;
		for (auto val : mse)
			sum += val;

		cout << endl;*/
		/*for (int i = 0; i < GROUP_SIZE; i++)
		{
			cout << mse[i]/variance << " ";
		}
		cout << endl;*/

	int min = 0, max = 0, sec_min = 0;
	int max_ = 0, min_ = 0;
	int max_num = 0, min_num = 1e10;
	double max_val = 0, min_val = 1e10, sec_min_val = 1e10;


	for (int i = 0; i < GROUP_SIZE; i++)
	{
		if (mse[i] < min_val)
		{
			min = i; min_val = mse[i];
		}
		if (mse[i] > max_val)
		{
			max = i; max_val = mse[i];
		}
		if (groupSize[i] < min_num)
		{
			min_ = i; min_num = groupSize[i];
		}
		if (groupSize[i] > max_num)
		{
			max_ = i; max_num = groupSize[i];
		}
	}
	for (int i = 0; i < GROUP_SIZE && i != min; i++)
	{
		if (mse[i] < sec_min_val)
		{
			sec_min = i; sec_min_val = mse[i];
		}
	}


	if (max_val > 4 * min_val)
	{
		//for (int i = 0; i < GROUP_SIZE; i++)
		//	cout << groupSize[i] << " ";

		//printinitial(rx, centroid);
		centroid[min][0] = centroid[max][0] + DELTA;
		centroid[min][1] = centroid[max][1] + DELTA;
		centroid[max][0] -= DELTA;
		centroid[max][1] -= DELTA;
		//printinitial(rx, centroid);
		Grouping(rx, centroid, group, groupSize);
		centroid[sec_min][0] = centroid[sec_min][1] = 0;
		centroid[min][0] = centroid[min][1] = 0;
		centroid[max][0] = centroid[max][1] = 0;


		for (int i = 0; i < groupSize[sec_min]; i++)
		{
			centroid[sec_min][0] += rx[group[sec_min][i]][0];
			centroid[sec_min][1] += rx[group[sec_min][i]][1];
		}
		centroid[sec_min][0] /= groupSize[sec_min];
		centroid[sec_min][1] /= groupSize[sec_min];
		for (int k = 0; k < groupSize[min]; k++)
		{
			centroid[min][0] += rx[group[min][k]][0];
			centroid[min][1] += rx[group[min][k]][1];
		}
		centroid[min][0] /= groupSize[min];
		centroid[min][1] /= groupSize[min];
		for (int k = 0; k < groupSize[max]; k++)
		{
			centroid[max][0] += rx[group[max][k]][0];
			centroid[max][1] += rx[group[max][k]][1];
		}
		centroid[max][0] /= groupSize[max];
		centroid[max][1] /= groupSize[max];


		//for (int i = 0; i < GROUP_SIZE; i++)
		//	cout << groupSize[i] << " ";

		return false;
	}

	return true;

}

void Inner_K_means(int min, int max, int sec_min, int* groupSize, int** group, double** centroid, double** rx)
{
	////////////////- inner grouping
	//cout << max << " " << min << " " << sec_min;
	//system("pause");

	double variation = 1e10;
	double** centroid_temp = new double* [2];
	int Cluster_Num = groupSize[max];
	for (int i = 0; i < 2; i++)
		centroid_temp[i] = new double[2];

	for (int iteration = 0;; iteration++)
	{
		//	printinitial(rx, centroid);
		vector<int> temp(Cluster_Num);
		for (int i = 0; i < groupSize[max]; i++)
			temp[i] = group[max][i];

		if (iteration != 0)
		{
			for (int i = 0; i < groupSize[min]; i++)
				temp[i + groupSize[max]] = group[min][i];
		}

		groupSize[min] = 0;
		groupSize[max] = 0;
		for (int i = 0; i < Cluster_Num; i++)
		{
			if (EuclideanDistance(rx[temp[i]], centroid[min]) < EuclideanDistance(rx[temp[i]], centroid[max]))
				group[min][groupSize[min]++] = temp[i];
			else
				group[max][groupSize[max]++] = temp[i];
		}

		//Grouping(rx, centroid, group, groupSize);
		//////////////////

		//////////////////-centroid renew

		centroid[min][0] = centroid[min][1] = 0;
		centroid[max][0] = centroid[max][1] = 0;

		for (int k = 0; k < groupSize[min]; k++)
		{
			centroid[min][0] += rx[group[min][k]][0];
			centroid[min][1] += rx[group[min][k]][1];
		}
		centroid[min][0] /= groupSize[min];
		centroid[min][1] /= groupSize[min];
		for (int k = 0; k < groupSize[max]; k++)
		{
			centroid[max][0] += rx[group[max][k]][0];
			centroid[max][1] += rx[group[max][k]][1];
		}
		centroid[max][0] /= groupSize[max];
		centroid[max][1] /= groupSize[max];



		variation = sqrt(EuclideanDistance(centroid_temp[0], centroid[min]) + EuclideanDistance(centroid_temp[1], centroid[max]));
		//cout << variation << " -  -";
		//printinitial(rx, centroid);
		if (iteration != 0 && variation < 1e-05)
			break;

		centroid_temp[0][0] = centroid[min][0]; centroid_temp[0][1] = centroid[min][1];
		centroid_temp[1][0] = centroid[max][0]; centroid_temp[1][1] = centroid[max][1];
		//////////////////////////////
	}

	Grouping(rx, centroid, group, groupSize);
	CentroidRenewal2(rx, centroid, group, groupSize);

	for (int i = 0; i < 2; i++)
		delete[] centroid_temp[i];
	delete[]centroid_temp;
}

bool ConditionCheck1(double** centroid, int* groupSize, double variance, double** rx, double** estimate, double** chCoef, int& count, int** group)
{
	/////////////////// -check symetric 
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
	if (sqrt(average) < 0.05 + count * 0.001) // ad-hoc
	{
		return false;
	}
	/////////////////////


	//printinitial(rx, centroid); 

	///////////////////// -modified method
	vector<double> mse(GROUP_SIZE);

	for (int i = 0; i < GROUP_SIZE; i++)
	{
		for (int j = 0; j < groupSize[i]; j++)
		{
			mse[i] += 0.5 * (pow(centroid[i][0] - rx[group[i][j]][0], 2) + pow(centroid[i][1] - rx[group[i][j]][1], 2));
		}
		mse[i] /= groupSize[i];
	}


	int min = 0, max = 0, sec_min = 0;
	int max_ = 0, min_ = 0;
	int max_num = 0, min_num = 1e10;
	double max_val = 0, min_val = 1e10, sec_min_val = 1e10;


	for (int i = 0; i < GROUP_SIZE; i++)
	{
		if (mse[i] < min_val)
		{
			min = i; min_val = mse[i];
		}
		if (mse[i] > max_val)
		{
			max = i; max_val = mse[i];
		}
		if (groupSize[i] < min_num)
		{
			min_ = i; min_num = groupSize[i];
		}
		if (groupSize[i] > max_num)
		{
			max_ = i; max_num = groupSize[i];
		}
	}
	for (int i = 0; i < GROUP_SIZE && i != min; i++)
	{
		if (mse[i] < sec_min_val)
		{
			sec_min = i; sec_min_val = mse[i];
		}
	}



	if (max_val > variance * min_val)
	{
		//printinitial(rx, centroid);
		centroid[min][0] = centroid[max][0] + DELTA;
		centroid[min][1] = centroid[max][1] + DELTA;
		centroid[max][0] -= DELTA;
		centroid[max][1] -= DELTA;
		//printinitial(rx, centroid);

		Inner_K_means(min, max, sec_min, groupSize, group, centroid, rx);

	}

	//printinitial(rx, centroid);
	// -check symetric 
	for (int i = 0; i < GROUP_SIZE; i++)
	{
		if (groupSize[i] == 0) return true; // empty cluster
	}
	average = 0;
	for (int i = 0; i < 2; i++)
	{
		double temp = 0;
		for (int j = 0; j < GROUP_SIZE; j++)
		{
			temp += centroid[j][i];
		}
		average += pow(temp / GROUP_SIZE, 2);
	}
	if (sqrt(average) < 0.05 + count * 0.001) // ad-hoc
	{
		return false;
	}
	count++;
	return true;
}

void CentroidRenewal2(double** rx, double** centroid, int** group, int* groupSize)
{
	for (int i = 0; i < GROUP_SIZE; i++)
	{
		for (int j = 0; j < 2; j++) // real and imaginary
		{
			double temp = 0;
			for (int k = 0; k < groupSize[i]; k++)
			{
				temp += rx[group[i][k]][j];
			}
			temp /= groupSize[i];							//�s��centroid
			centroid[i][j] = temp;
		}
		//	cout << variation[i] << " ";
	}
	//	cout << endl;
}