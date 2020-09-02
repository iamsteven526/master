#define _USE_MATH_DEFINES
#include <iostream>
#include <random>
#include "parameters.h"
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <cstring>
#include <cmath>
#pragma warning(disable:4996)

using namespace std;

void Clustering(double **rx, double **centroid, int **group, int *groupSize, double *distList, double *variation, double **softAssign, double variance, double **estimate, double ***chCoef)
{
	bool reset = true;
	int count = 0;

	//---- generater paired samples
	for (int i = 0; i < CODE_LEN; i++)
	{
		rx[i + CODE_LEN][0] = -rx[i][0];
		rx[i + CODE_LEN][1] = -rx[i][1];
	}


	while (true)
	{
		count++;
		if (reset)
		{
			if (PROPOSAL == 0)
				InitialSeeding(rx, centroid, group, groupSize, distList, variance, NUM_USER, GROUP_SIZE, CODE_LEN);
			else
				InitialSeeding(rx, centroid, group, groupSize, distList, variance, NUM_USER, GROUP_SIZE, 2 * CODE_LEN);
			reset = false;
		}

		if (PROPOSAL == 0)
			Grouping(rx, centroid, group, groupSize, NUM_USER, GROUP_SIZE, CODE_LEN);
		else
			Grouping(rx, centroid, group, groupSize, NUM_USER, GROUP_SIZE, 2 * CODE_LEN);
		CentroidRenewal(rx, centroid, group, groupSize, variation, GROUP_SIZE);

		bool convergence = true;
		for (int i = 0; i < GROUP_SIZE; i++)
		{
			if (variation[i] > 1e-5)
			{
				convergence = false;
				break;
			}
		}


		if (convergence == 1)
		{
			if (PROPOSAL == 0) //---- Original
			{
				reset = ConditionCheck(centroid, groupSize, count, GROUP_SIZE);
				//reset = ConditionCheck1(centroid, groupSize, variance, rx, estimate, chCoef, count, group);

				if (!reset)
				{
					if (EM_GMM) EMClustering(rx, centroid, softAssign, variance, CODE_LEN, GROUP_SIZE);
					CoefEstimation(centroid, estimate, variance, reset);
					return;
				}
			}
			else if (PROPOSAL == 1) //---- Proposal-1
			{
				int inner_count = 0; //---- count for inner k-means

				for (;;)
				{
					//cout << count << " ";
					reset = 0;
					//---- symmetric check
					for (int i = 0; i < GROUP_SIZE; i++)
					{
						if (groupSize[i] == 0)
						{
							reset = true; // empty cluster
						}
					}
					if (reset)
						break;

					double average = 0;
					for (int i = 0; i < 2; i++)
					{
						double temp = 0;
						for (int j = 0; j < GROUP_SIZE; j++)
						{
							temp += centroid[j][i];
						}
						average += pow(temp, 2);
					}
					//---- 
					if (0.5 * average > variance + count * 0.01)
						reset = 1;

					if (reset)
					{
						//---- mse calculating
						vector<double> mse(GROUP_SIZE);
						for (int i = 0; i < GROUP_SIZE; i++)
						{
							for (int j = 0; j < groupSize[i]; j++)
							{
								mse[i] += 0.5 * (pow(centroid[i][0] - rx[group[i][j]][0], 2) + pow(centroid[i][1] - rx[group[i][j]][1], 2));
							}
							mse[i] /= groupSize[i];
						}

						//---- Sorting group size
						int max_tr = 0, min_tr = 0, sec_min_tr = 0, max_val = 0, min_val = 1e10, sec_min_val = 1e10;

						for (int i = 0; i < GROUP_SIZE; i++)
						{
							if (mse[i] < min_val)
							{
								min_tr = i; min_val = mse[i];
							}
							if (mse[i] > max_val)
							{
								max_tr = i; max_val = mse[i];
							}
						}
						for (int i = 0; i < GROUP_SIZE && i != min_tr; i++)
						{
							if (mse[i] < sec_min_val)
							{
								sec_min_tr = i; sec_min_val = mse[i];
							}
						}

						//---- centroid delte and centroid split
						if (max_val > variance * min_val)
						{
							centroid[min_tr][0] = centroid[max_tr][0] + DELTA;
							centroid[min_tr][1] = centroid[max_tr][1] + DELTA;
							centroid[max_tr][0] -= DELTA;
							centroid[max_tr][1] -= DELTA;

							Inner_K_means(min_tr, max_tr, sec_min_tr, groupSize, group, centroid, rx);
						}
						else
							break;
					}
					else
						break;

					inner_count++;

					//---- if inner_count > 100 reset for the new initial centroids
					if (inner_count > 100 )
					{
						break;
					}

				}

				if (count > 1000)
				{

					int inner_count = 0;  //---- count for inner k-means

					for (;;)
					{
						int max_tr = 0, min_tr = 0, sec_min_tr = 0, max_val = 0, min_val = 1e10, sec_min_val = 1e10;

						//---- Sorting group size
						for (int i = 0; i < GROUP_SIZE; i++)
						{
							if (groupSize[i] < min_val)
							{
								min_tr = i; min_val = groupSize[i];
							}
							if (groupSize[i] > max_val)
							{
								max_tr = i; max_val = groupSize[i];
							}
						}
						for (int i = 0; i < GROUP_SIZE && i != min_tr; i++)
						{
							if (groupSize[i] < sec_min_val)
							{
								sec_min_tr = i; sec_min_val = groupSize[i];
							}
						}

						//---- centroid delte and centroid split
						if (max_val > 2 * min_val)
						{
							centroid[min_tr][0] = centroid[max_tr][0] + DELTA;
							centroid[min_tr][1] = centroid[max_tr][1] + DELTA;
							centroid[max_tr][0] -= DELTA;
							centroid[max_tr][1] -= DELTA;

							Inner_K_means(min_tr, max_tr, sec_min_tr, groupSize, group, centroid, rx);
						}
						else
							break;
						inner_count++;

						//---- if inner_count > 100 reset for the new initial centroids
						if (inner_count > 100)
						{
							count++;
							break;
						}
					}
					if (EM_GMM) EMClustering(rx, centroid, softAssign, variance, 2 * CODE_LEN, GROUP_SIZE);
					CoefEstimation(centroid, estimate, variance, reset);
					return;
				}

				if (!reset)
				{
					if (EM_GMM) EMClustering(rx, centroid, softAssign, variance, 2 * CODE_LEN, GROUP_SIZE);
					CoefEstimation(centroid, estimate, variance, reset);

					if (!reset)
						return;
				}
			}
			else //---- Proposal-2
			{
				int inner_count = 0;  //---- count for inner k-means
				reset = 0;
				for (;;)
				{
					int max_tr = 0, min_tr = 0, sec_min_tr = 0, max_val = 0, min_val = 1e10, sec_min_val = 1e10;

					//---- Sorting group size
					for (int i = 0; i < GROUP_SIZE; i++)
					{
						if (groupSize[i] < min_val)
						{
							min_tr = i; min_val = groupSize[i];
						}
						if (groupSize[i] > max_val)
						{
							max_tr = i; max_val = groupSize[i];
						}
					}
					for (int i = 0; i < GROUP_SIZE && i != min_tr; i++)
					{
						if (groupSize[i] < sec_min_val)
						{
							sec_min_tr = i; sec_min_val = groupSize[i];
						}
					}

					//---- centroid delte and centroid split
					if (max_val > 2 * min_val)
					{
						centroid[min_tr][0] = centroid[max_tr][0] + DELTA;
						centroid[min_tr][1] = centroid[max_tr][1] + DELTA;
						centroid[max_tr][0] -= DELTA;
						centroid[max_tr][1] -= DELTA;

						Inner_K_means(min_tr, max_tr, sec_min_tr, groupSize, group, centroid, rx);
					}
					else
						break;
					inner_count++;

					//---- if inner_count > 100 reset for the new initial centroids
					if (inner_count > 100)
					{
						count++;
						break;
					}
				}

				//---- count is the reset times 
				if (!reset || count > 100)
				{
					if (EM_GMM) EMClustering(rx, centroid, softAssign, variance, 2 * CODE_LEN, GROUP_SIZE);
					CoefEstimation(centroid, estimate, variance, reset);

					return;
				}

			}
		}
	}
}

void InitialSeeding(double** rx, double** centroid, int** group, int* groupSize, double* distList, double variance, int CLUSTER_USER, int CLUSTER_GROUP, int CLUSTER_LEN)
{
	random_device seed;
	mt19937 generator(seed());
	uniform_real_distribution<double> uniform(0, 1);
	if (INI_METHOD == 0) // LBG
	{
		double delta = 0.001;
		double phi = 2 * M_PI * uniform(generator); // splitting phase
		vector<double> dist(CLUSTER_GROUP);
		centroid[0][0] = +delta * cos(phi); centroid[0][1] = +delta * sin(phi);
		centroid[1][0] = -delta * cos(phi); centroid[1][1] = -delta * sin(phi);
		for (int i = 1; i < CLUSTER_USER; i++)
		{
			memset(groupSize, 0, sizeof(int) * CLUSTER_GROUP);
			for (int j = 0; j < CLUSTER_LEN; j++) // grouping
			{
				for (int k = 0; k < (1 << i); k++)
				{
					dist[k] = EuclideanDistance(rx[j], centroid[k]);
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
					centroid[2 * j][0] += rx[group[j][k]][0];
					centroid[2 * j][1] += rx[group[j][k]][1];
				}
				centroid[2 * j][0] /= groupSize[j]; centroid[2 * j + 1][0] = centroid[2 * j][0];
				centroid[2 * j][1] /= groupSize[j]; centroid[2 * j + 1][1] = centroid[2 * j][1];
				centroid[2 * j + 0][0] += delta * cos(phi); centroid[2 * j + 0][1] += delta * sin(phi);
				centroid[2 * j + 1][0] -= delta * cos(phi); centroid[2 * j + 1][1] -= delta * sin(phi);
			}
		}

	}
	else if (INI_METHOD == 1 || INI_METHOD == 2) // k-means++
	{
		//uniform_int_distribution<int> ranPick(0, CLUSTER_LEN - 1);
		//int iniPick = ranPick(generator);
		int iniPick = rand() % CLUSTER_LEN;
		centroid[0][0] = rx[iniPick][0]; centroid[0][1] = rx[iniPick][1];
		for (int i = 1; i < CLUSTER_GROUP; i++)
		{
			double sum = 0;
			for (int j = 0; j < CLUSTER_LEN; j++) // distance from the closest centroid
			{
				distList[j] = EuclideanDistance(rx[j], centroid[0]);
				for (int k = 1; k < i; k++)
				{
					double dist = EuclideanDistance(rx[j], centroid[k]);
					if (dist < distList[j]) distList[j] = dist;
				}
				distList[j] = pow(distList[j], 2);
				sum += distList[j];
			}
			if (INI_METHOD == 2)
			{
				vector<int> sort_seq(CLUSTER_LEN);
				vector<double> sort_val(CLUSTER_LEN);
				int s = 0;
				generate(sort_seq.begin(), sort_seq.end(), [&] {return s++; });
				for (int l = 0; l < CLUSTER_LEN; l++)
				{
					sort_val[l] = distList[l];
				}
				sort(sort_seq.begin(), sort_seq.end(), [&](int x, int y) {return sort_val[x] < sort_val[y]; });


				for (int j = 0; j < CLUSTER_LEN - (CLUSTER_LEN / CLUSTER_GROUP); j++)
				{
					sum -= distList[sort_seq[j]];
					distList[sort_seq[j]] = 0;
				}
			}
			double ranNum = uniform(generator) * sum;
			for (int j = 0; j < CLUSTER_LEN; j++) // select with weighted probability
			{
				ranNum -= distList[j];
				if (ranNum <= 0)
				{
					centroid[i][0] = rx[j][0]; centroid[i][1] = rx[j][1];
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

void Grouping(double** rx, double** centroid, int** group, int* groupSize, int CLUSTER_USER, int CLUSTER_GROUP, int CLUSTER_LEN)
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

	/*
	for (int i = 0; i < groupSize[0]; i++)
	{
		cout << rx[group[0][i]][0]<<" ";
	}
	cout << endl;
	for (int i = 0; i < groupSize[0]; i++)
	{
		cout << rx[group[0][i]][1] << " ";
	}
	cout << endl;
	*/

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

	Grouping(rx, centroid, group, groupSize, NUM_USER, GROUP_SIZE, CODE_LEN);
	CentroidRenewal2(rx, centroid, group, groupSize);
}

double EuclideanDistance(double *x, double *y)
{
	return sqrt(pow(x[0] - y[0], 2) + pow(x[1] - y[1], 2));
}

void CentroidRenewal(double** rx, double** centroid, int** group, int* groupSize, double* variation, int CLUSTER_GROUP)
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
			temp /= groupSize[i];							//�s��centroid
			variation[i] += pow(centroid[i][j] - temp, 2);	//�s�¤������t�O
			centroid[i][j] = temp;
		}
		variation[i] = sqrt(variation[i]);
		//	cout << variation[i] << " ";
	}
	//	cout << endl;
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

void CoefEstimation(double** centroid, double** estimate, double variance, bool& reset)
{
	int num_level = NUM_LEVEL;
	int num_user = NUM_USER;
	vector<vector<double>> sup_centroid(num_level, vector<double>(2));
	for (int i = 0; i < num_level; i++)
		for (int j = 0; j < 2; j++)
			sup_centroid[i][j] = centroid[i][j];

	// Until all chCoefs are derived
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{

		vector<vector<double>> pair(num_level / 2, vector<double>(2)); // pair and inverse pair 
		vector<double> group_centroid(2);								// group centroid of pair	
		vector<double> inverse_group_centroid(2);						// group centroid of inverse pair
		vector<int> pair_num(num_level);								// the number of centroids which need to pair, first step seq
		vector<int> pair_num_;											// the number of centroids which need to pair, sec step seq
		vector<int> inverse_pair_num_;									// the inverse of pair_num_

		int p = 0;
		generate(pair_num.begin(), pair_num.end(), [&] {return p++; });
		double mse = 0;

		//---find pair

		pair = pair_seq(sup_centroid, num_level, mse, group_centroid, pair_num);

		if (variance < mse)
			reset = 1;


		if (reset)
			break;

		if (num_level != 2)
		{
			//---find channel coefficient
			for (int i = 0; i < pow(2, num_level / 2); i++)
			{
				//--- list all state of pair arrangement
				int reg = i;
				for (int j = 0; j < num_level / 2; j++)
				{
					int num = pair[j][reg % 2];
					group_centroid[0] += sup_centroid[num][0];
					group_centroid[1] += sup_centroid[num][1];
					pair_num_.push_back(pair[j][reg % 2]);
					inverse_pair_num_.push_back(pair[j][(reg + 1) % 2]);
					inverse_group_centroid[0] += sup_centroid[pair[j][(reg + 1) % 2]][0];
					inverse_group_centroid[1] += sup_centroid[pair[j][(reg + 1) % 2]][1];

					reg = reg / 2;
				}
				group_centroid[0] /= num_level / 2; group_centroid[1] /= num_level / 2;
				inverse_group_centroid[0] /= num_level / 2; inverse_group_centroid[1] /= num_level / 2;

				//---check symmtric
				vector<vector<double>> pair_(num_level / 4);
				pair_ = pair_seq(sup_centroid, num_level / 2, mse, group_centroid, pair_num_);

				//--- Derive correct chCoef
				if (variance > mse)
				{
					vector<vector<double>> temp(num_level / 2, vector<double>(2));

					for (int j = 0; j < num_level / 2; j++)
					{
						temp[j][0] = (sup_centroid[pair_num_[j]][0] - group_centroid[0] - (sup_centroid[inverse_pair_num_[j]][0] - inverse_group_centroid[0])) / 2;
						temp[j][1] = (sup_centroid[pair_num_[j]][1] - group_centroid[1] - (sup_centroid[inverse_pair_num_[j]][1] - inverse_group_centroid[1])) / 2;
					}

					for (int j = 0; j < num_level / 2; j++)
						for (int k = 0; k < 2; k++)
							sup_centroid[j][k] = temp[j][k];

					estimate[nuser][0] = (group_centroid[0] - inverse_group_centroid[0]) / 2;
					estimate[nuser][1] = (group_centroid[1] - inverse_group_centroid[1]) / 2;
					num_level /= 2;
					break;
				}

				if (variance < mse && i + 1 == pow(2, num_level / 2))
				{
					//cout << "*";
					reset = 1;
					break;
					//	cout << "camn't find symmetric centroid";
					//	system("pause");


				}
				//---clear previous state
				pair_num_.clear();
				inverse_pair_num_.clear();
				group_centroid[0] = 0; group_centroid[1] = 0;
				inverse_group_centroid[0] = 0; inverse_group_centroid[1] = 0;

				//system("pause");
			}
		}
		else
		{
			estimate[nuser][0] = (sup_centroid[0][0] - sup_centroid[1][0]) / 2;
			estimate[nuser][1] = (sup_centroid[0][1] - sup_centroid[1][1]) / 2;
		}
	}
}

void MSEComparison(double*** chCoef, double** estimate, double& mse)
{
	unordered_set<int> map;
	double temp[NUM_USER][2] = { 0 };

	double min_value;
	int min_tr;


	for (int i = 0; i < NUM_USER; i++)
	{
		min_value = 1000;
		min_tr = 0;
		for (int j = 1; j <= NUM_USER; j++)
		{
			if (map.count(j))
				continue;
			double reg1 = pow(chCoef[i][0][0] - estimate[j - 1][0], 2) + pow(chCoef[i][1][0] - estimate[j - 1][1], 2);
			double reg2 = pow(chCoef[i][0][0] + estimate[j - 1][0], 2) + pow(chCoef[i][1][0] + estimate[j - 1][1], 2);

			if (reg1 < reg2)
			{
				if (reg1 < min_value)
				{
					min_value = reg1;
					min_tr = j;
				}
			}
			else
			{
				if (reg2 < min_value)
				{
					min_value = reg2;
					min_tr = -j;
				}
			}
		}

		map.insert(min_tr > 0 ? min_tr : -min_tr);
		temp[i][0] = min_tr > 0 ? estimate[min_tr - 1][0] : -estimate[-min_tr - 1][0];
		temp[i][1] = min_tr > 0 ? estimate[min_tr - 1][1] : -estimate[-min_tr - 1][1];
		//mse += 0.5 * min_value;
	}

	for (int i = 0; i < NUM_USER; i++)
	{
		estimate[i][0] = temp[i][0];
		estimate[i][1] = temp[i][1];
	}

	
}

bool ConditionCheck(double** centroid, int* groupSize, int& count, int CLUSTER_GROUP)
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
	if (sqrt(average) > 0.05 + count * 0.001) // ad-hoc
	{
		//count++;
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

vector<vector<double>> pair_seq(vector<vector<double>> centroid, int num_level, double &mse, vector<double> group_centroid, vector<int> pair_num)
{
	unordered_set<int> map;
	vector<vector<double>> pair(num_level/2,vector<double> (2));

	mse = 0;
	int m = 0;
	for (int i = 0; i < num_level; i++)
	{
		if (map.count(pair_num[i]) == 1)
			continue;
		map.insert(pair_num[i]);
		double min = 1000;
		int min_ = 0;
		for (int j = 0; j < num_level; j++)
		{
			if (map.count(pair_num[j]) == 1)
				continue;
			double reg = pow(centroid[pair_num[i]][0] + centroid[pair_num[j]][0] - 2 * group_centroid[0], 2) + pow(centroid[pair_num[i]][1] + centroid[pair_num[j]][1] - 2 * group_centroid[1], 2);
			if (reg < min)
			{
				min = reg;
				min_ = pair_num[j];
			}
		}
		mse += min*0.5;
		map.insert(min_);
		pair[m][0] = pair_num[i];
		pair[m][1] = min_;
		m++;
	}
	
	mse /= num_level / 2;
	return pair;
}

void printinitial(double** rx, double** centroid, double*** chCoef)
{
	FILE* rx_x_intitial = fopen("rx_x_intitial.txt", "w");
	FILE* rx_y_intitial = fopen("rx_y_intitial.txt", "w");
	FILE* centroid_x_intitial = fopen("centroid_x_intitial.txt", "w");
	FILE* centroid_y_intitial = fopen("centroid_y_intitial.txt", "w");
	FILE* chCoef_x_intitial = fopen("chCoef_x_intitial.txt", "w");
	FILE* chCoef_y_intitial = fopen("chCoef_y_intitial.txt", "w");
	for (int j = 0; j < CODE_LEN; j++)
	{
		//cout << rx[j][0] << " ";
		fprintf(rx_x_intitial, "%e ", rx[j][0]);
		
	}

	for (int j = 0; j < CODE_LEN; j++)
	{
		fprintf(rx_y_intitial, "%e ", rx[j][1]);
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

	for (int j = 0; j < NUM_USER; j++)
	{
		fprintf(chCoef_x_intitial, "%e ", chCoef[j][0][0]);
		fprintf(chCoef_x_intitial, "%e ", -chCoef[j][0][0]);
	}

	for (int j = 0; j < NUM_USER; j++)
	{
		fprintf(chCoef_y_intitial, "%e ", chCoef[j][1][0]);
		fprintf(chCoef_y_intitial, "%e ", -chCoef[j][1][0]);
	}

	fclose(rx_x_intitial);
	fclose(rx_y_intitial);
	fclose(centroid_x_intitial);
	fclose(centroid_y_intitial);
	fclose(chCoef_x_intitial);
	fclose(chCoef_y_intitial);

	system("pause");
}

