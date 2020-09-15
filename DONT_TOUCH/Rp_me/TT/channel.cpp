#include<iostream>
#include<random>
#include "parameter.h"

using namespace std;

random_device seed;
mt19937 generator(seed());
normal_distribution<double> normal(0, 1);

void channel(double stdDev, vector<int>&tx, vector<double>&rx,int recursive)
{
	if (recursive < 1)
	{
		for (int i = 0; i < LDPC_H_COL; i++)
		{
			rx[recursive*IR_LEN + i] = 1 - (2 * tx[i]) + stdDev * normal(generator);
		//	rx[recursive*IR_LEN + i] = 1 - (2 * tx[i]);
		}
	}
	else 
	{
		for (int i = 0; i < IR_LEN; i++)
		{
			rx[LDPC_H_COL +(recursive-1)*IR_LEN + i] = 1 - (2 * tx[i]) + stdDev * normal(generator);
		//	rx[LDPC_H_COL + (recursive - 1)*IR_LEN + i] = 1 - (2 * tx[i]);
		}
	}
}