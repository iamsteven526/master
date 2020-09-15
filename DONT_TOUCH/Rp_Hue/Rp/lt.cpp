#include <iostream>
#include <random>
#include "lt.h"
#include "parameters.h"
#include "degree_profile.h"
using namespace std;

extern random_device seed;

LT::LT(int input_len, int max_output_len, int incremental_len)
{
	mt19937 generator(seed());
	uniform_int_distribution<int> uniform(0, RAND_MAX);
	in_len = input_len;
	max_out_len = max_output_len;
	ir_len = incremental_len;
	o_deg = new int[max_out_len];			// degree of output check nodes 
	i_deg = new int[in_len];				// degree of input nodes
	o_edge = new int*[max_out_len];			// connected input nodes of each output check node
	i_edge = new int*[in_len];				// connected output check nodes of each input node
	o_order = new int*[max_out_len];		// order of message passing
	i_order = new int*[in_len];				// order of message passing
	llr_o2i = new double*[max_out_len];		// LLRs from output check nodes to input nodes
	llr_i2o = new double*[in_len];			// LLRs from input nodes to output check nodes
	ch_llr = new double[max_out_len];		// channel LLRs
	ch_value = new double[max_out_len];		// tanh( ch_value / 2 )
	//---------- calculate the cdf of output check node degree distribution ----------
	double **deg_cdf = new double *[STEP_SIZE];
	for (int i = 0; i < STEP_SIZE; i++)
	{
		deg_cdf[i] = new double[DEG_NUM];
	}
	for (int i = 0; i < STEP_SIZE; i++)
	{
		deg_cdf[i][0] = 0;
		for (int j = 1; j < DEG_NUM; j++)
		{
			deg_cdf[i][j] = deg_cdf[i][j - 1] + deg_pdf[i][j];
		}
		for (int j = 1; j < DEG_NUM; j++)
		{
			deg_cdf[i][j] = deg_cdf[i][j] / deg_cdf[i][DEG_NUM - 1]; // normalization
		} 
	}
	//---------- set up the degree and edge for both input node and output check node ----------
	memset(i_deg, 0, in_len * sizeof(int));
	if (SYSTEMATIC)
	{
		//---------- select the degree of output check nodes ----------
		int step_count = 0;
		for (int i = 0; i < in_len; i++)
		{
			o_deg[i] = 1; // systematic
			i_deg[i]++;
		}
		for (int i = in_len; i < max_out_len; i++)
		{
			if ((i >= step[step_count]) && (STEP_SIZE > 1) && (step_count < STEP_SIZE - 1)) step_count++; // switch to next degree distribution
			double random = (double)(uniform(generator) + 1) / (double)(RAND_MAX + 2);
			for (int j = 1; j < DEG_NUM; j++)
			{
				if (deg_cdf[step_count][j - 1] < random && random <= deg_cdf[step_count][j])
				{
					o_deg[i] = deg[j];
					break;
				}
			}
		}
		//---------- set up the edge connection from output check nodes ----------
		for (int i = 0; i < max_out_len; i++)
		{
			o_edge[i] = new int[o_deg[i]];
			o_order[i] = new int[o_deg[i]];
			llr_o2i[i] = new double[o_deg[i]];
		}
		for (int i = 0; i < max_out_len; i++)
		{
			if (i < in_len) o_edge[i][0] = i; // systematic 
			else
			{
				for (int j = 0; j < o_deg[i]; j++)
				{
					bool flag = true;
					while (flag)
					{
						o_edge[i][j] = uniform(generator) % in_len;
						flag = false;
						for (int k = 0; k < j; k++)
						{
							if (o_edge[i][k] == o_edge[i][j])
							{
								flag = true;
								break; // select again
							}
						}
					}
					i_deg[o_edge[i][j]]++;
				}
			}
		}
	}
	else // non-systematic
	{
		//---------- select the degree of output check nodes ----------
		int step_count = 0;
		for (int i = 0; i < max_out_len; i++)
		{
			if ((i >= step[step_count]) && (STEP_SIZE > 1) && (step_count < STEP_SIZE - 1)) step_count++; // switch to next degree distribution
			double random = (double)(uniform(generator) + 1) / (double)(RAND_MAX + 2);
			for (int j = 1; j < DEG_NUM; j++)
			{
				if (deg_cdf[step_count][j - 1] < random && random <= deg_cdf[step_count][j])
				{
					o_deg[i] = deg[j];
					break;
				}
			}
		}
		//---------- set up the edge connection from output check nodes ----------
		for (int i = 0; i < max_out_len; i++)
		{
			o_edge[i] = new int[o_deg[i]];
			o_order[i] = new int[o_deg[i]];
			llr_o2i[i] = new double[o_deg[i]];
		}
		for (int i = 0; i < max_out_len; i++)
		{
			for (int j = 0; j < o_deg[i]; j++)
			{
				bool flag = true;
				while (flag)
				{
					o_edge[i][j] = uniform(generator) % in_len;
					flag = false;
					for (int k = 0; k < j; k++)
					{
						if (o_edge[i][k] == o_edge[i][j])
						{
							flag = true;
							break; // select again
						}
					}
				}
				i_deg[o_edge[i][j]]++;
			}
		}
	}
	//---------- set up the edge connection from input nodes ----------
	int *i_count = new int[in_len];
	for (int i = 0; i < in_len; i++)
	{
		i_edge[i] = new int[i_deg[i]];
		i_order[i] = new int[i_deg[i]];
		llr_i2o[i] = new double[i_deg[i]];
	}
	memset(i_count, 0, in_len * sizeof(int));
	for (int i = 0; i < max_out_len; i++)
	{
		for (int j = 0; j < o_deg[i]; j++)
		{
			i_edge[o_edge[i][j]][i_count[o_edge[i][j]]] = i;
			i_order[o_edge[i][j]][i_count[o_edge[i][j]]] = j;
			o_order[i][j] = i_count[o_edge[i][j]];
			i_count[o_edge[i][j]]++;
		}
	}
	for (int i = 0; i < STEP_SIZE; i++)
	{
		delete[] deg_cdf[i];
	}
	delete[] deg_cdf;
	delete[] i_count;
}

void LT::Encoder(int *data, int *tx)
{
	for (int i = 0; i < max_out_len; i++)
	{
		tx[i] = 0;
		for (int j = 0; j < o_deg[i]; j++)
		{
			tx[i] ^= data[o_edge[i][j]];
		}
	}
}

void LT::Decoder(double variance, double *rx, double *app_llr, int current_len, int it_max)
{
	

	for (int i = 0; i < current_len; i++)
	{
		ch_llr[i] = (2 / variance)*rx[i];
		if (ch_llr[i] >= 30) { ch_value[i] = 1;}
		else if (ch_llr[i] <= -30) ch_value[i] = -1;
		else ch_value[i] = tanh(ch_llr[i] / 2);
	}
	for (int i = 0; i < current_len; i++)
	{
		memset(llr_o2i[i], 0, sizeof(double)*o_deg[i]);
	}
	for (int i = 0; i < in_len; i++)
	{
		memset(llr_i2o[i], 0, sizeof(double)*i_deg[i]);
	}
	//---------- Belief propagation ----------
	for (int it = 1; it <= it_max; it++)
	{
		LLR_O2I(current_len);
		LLR_I2O(current_len);
	}
	memset(app_llr, 0, sizeof(double)*in_len);
	for (int i = 0; i < current_len; i++)
	{
		for (int j = 0; j < o_deg[i]; j++)
		{
			app_llr[o_edge[i][j]] += llr_o2i[i][j];
		}
	}

}

void LT::LLR_O2I(int current_len)
{
	for (int i = 0; i < current_len; i++)
	{
		if (o_deg[i] == 1) llr_o2i[i][0] = ch_llr[i]; 
		else
		{
			for (int j = 0; j < o_deg[i]; j++)
			{
				if (llr_i2o[o_edge[i][j]][o_order[i][j]] >= 30) { temp_llr[j] = 1;}
				else if (llr_i2o[o_edge[i][j]][o_order[i][j]] <= -30) temp_llr[j] = -1;
				else temp_llr[j] = tanh(llr_i2o[o_edge[i][j]][o_order[i][j]] / 2);
			}
			for (int j = 0; j < o_deg[i]; j++)
			{
				double temp = ch_value[i];

				for (int k = 0; k < o_deg[i]; k++)
				{
					if (k != j) temp *= temp_llr[k];
				}
				if (temp == 1) { llr_o2i[i][j] = 40;}
				else if (temp == 0) llr_o2i[i][j] = 0;
				else if (temp == -1) llr_o2i[i][j] = -40;
				else llr_o2i[i][j] = log((1 + temp) / (1 - temp));
			}
		}
	}
}

void LT::LLR_I2O(int current_len)
{
	for (int i = 0; i < in_len; i++)
	{
		for (int j = 0; j < i_deg[i]; j++)
		{
			if (i_edge[i][j] < current_len)
			{
				double temp = 0;
				for (int k = 0; k < i_deg[i]; k++)
				{
					if ((k != j) && (i_edge[i][k] < current_len)) temp += llr_o2i[i_edge[i][k]][i_order[i][k]];
					if (i_edge[i][k] >= current_len) break;
				}
				llr_i2o[i][j] = temp;
				
			}
			else break; 
		}
	}
}