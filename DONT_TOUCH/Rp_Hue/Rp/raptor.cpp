#include <iostream>
#include "ldpc.h"
#include "lt.h"
#include "parameters.h"
using namespace std;

void RaptorEncoder(LDPC &ldpc, LT &lt, int *data, int *inter_code, int *tx)
{
	for (int i = 0; i < DATA_LEN; i++)
	{
		data[i] = rand() % 2;
	}
	ldpc.Encoder(data, inter_code);
	lt.Encoder(inter_code, tx);
}

void TandemDecoder(LT &lt, LDPC &ldpc, double variance, double *rx, double *inter_llr, double *decoded_llr, int *data, int &ack, int &ir_count)
{
	for (ir_count = 0; lt.in_len + ir_count * lt.ir_len <= lt.max_out_len; ir_count++)
	{
		int current_len = lt.in_len + ir_count * lt.ir_len;
		lt.Decoder(variance, rx, inter_llr, current_len, LT_IT);
		ldpc.Decoder(variance, inter_llr, decoded_llr, LDPC_IT);
		ack = 1;
		for (int j = 0; j < ldpc.data_len; j++)
		{
			if (HARD(decoded_llr[j]) != data[j])
			{
				ack = 0; // retransmission
				break;
			}
		}
		if (ack == 1) return;
	}
	ir_count--;
}

void JointDecoder(LT &lt, LDPC &ldpc, double variance, double *rx, double *decoded_llr, int *data, int &ack, int &ir_count)
{
	for (ir_count = 0; lt.in_len + ir_count * lt.ir_len <= lt.max_out_len; ir_count++)
	{
		//---------- initialization ----------
		int current_len = lt.in_len + ir_count * lt.ir_len;
		for (int i = 0; i < ldpc.h_row; i++)
		{
			memset(ldpc.app_llr[i], 0, ldpc.h_col * sizeof(double));
		}
		for (int i = 0; i < current_len; i++)
		{
			lt.ch_llr[i] = (2 / variance)*rx[i];
			if (lt.ch_llr[i] >= 30) lt.ch_value[i] = 1;
			else if (lt.ch_llr[i] <= -30) lt.ch_value[i] = -1;
			else lt.ch_value[i] = tanh(lt.ch_llr[i] / 2);
		}
		for (int i = 0; i < current_len; i++)
		{
			memset(lt.llr_o2i[i], 0, sizeof(double)*lt.o_deg[i]);
		}
		for (int i = 0; i < lt.in_len; i++)
		{
			memset(lt.llr_i2o[i], 0, sizeof(double)*lt.i_deg[i]);
		}
		//---------- belief propagation ----------
		for (int it = 0; it < JOINT_IT; it++)
		{
			LLR_C2V(ldpc, lt, current_len); // check messages
			for (int i = 0; i < lt.in_len; i++)
			{
				decoded_llr[i] = 0;
				for (int j = 0; j < lt.i_deg[i]; j++)
				{
					if (lt.i_edge[i][j] < current_len) decoded_llr[i] += lt.llr_o2i[lt.i_edge[i][j]][lt.i_order[i][j]];
				}
				for (int j = 0; j < ldpc.col_w[i]; j++)
				{
					decoded_llr[i] += ldpc.ext_llr[ldpc.col_edge[i][j]][i];
				}
			}
			//---------- test ----------
			bool flag = false;
			for (int i = 0; i < ldpc.h_row; i++)
			{
				int temp = 0;
				for (int j = 0; j < ldpc.row_w[i]; j++)
				{
					temp ^= HARD(decoded_llr[ldpc.row_edge[i][j]]);
				}
				if (temp == 1)
				{
					flag = true;
					break;
				}
			}
			if (!flag) break; // early termination
			LLR_V2C(ldpc, lt, current_len); // bit messages
		}
		ack = 1;
		for (int j = 0; j < ldpc.data_len; j++)
		{
			if (HARD(decoded_llr[j]) != data[j])
			{
				ack = 0; // retransmission
				break;
			}
		}
		if (ack == 1) return;
	}
	ir_count--;
}

void LLR_C2V(LDPC &ldpc, LT &lt, int current_len) 
{
	for (int i = 0; i < current_len; i++)
	{
		if (lt.o_deg[i] == 1) lt.llr_o2i[i][0] = lt.ch_llr[i];
		else
		{
			for (int j = 0; j < lt.o_deg[i]; j++)
			{
				if (lt.llr_i2o[lt.o_edge[i][j]][lt.o_order[i][j]] >= 30) lt.temp_llr[j] = 1;
				else if (lt.llr_i2o[lt.o_edge[i][j]][lt.o_order[i][j]] == 0) lt.temp_llr[j] = 0;
				else if (lt.llr_i2o[lt.o_edge[i][j]][lt.o_order[i][j]] <= -30) lt.temp_llr[j] = -1;
				else lt.temp_llr[j] = tanh(lt.llr_i2o[lt.o_edge[i][j]][lt.o_order[i][j]] / 2);
			}
			for (int j = 0; j < lt.o_deg[i]; j++)
			{
				double temp = lt.ch_value[i];
				for (int k = 0; k < lt.o_deg[i]; k++)
				{
					if (k != j) temp *= lt.temp_llr[k];
				}
				if (temp == 1) lt.llr_o2i[i][j] = 40;
				else if (temp == 0) lt.llr_o2i[i][j] = 0;
				else if (temp == -1) lt.llr_o2i[i][j] = -40;
				else lt.llr_o2i[i][j] = log((1 + temp) / (1 - temp));
			}
		}
	}
	for (int j = 0; j < ldpc.h_row; j++)
	{
		for (int i = 0; i < ldpc.row_w[j]; i++)
		{
			if (ldpc.app_llr[j][ldpc.row_edge[j][i]] >= 30) ldpc.temp_llr[i] = 1;
			else if (ldpc.app_llr[j][ldpc.row_edge[j][i]] == 0) ldpc.temp_llr[i] = 0;
			else if (ldpc.app_llr[j][ldpc.row_edge[j][i]] <= -30) ldpc.temp_llr[i] = -1;
			else ldpc.temp_llr[i] = tanh(ldpc.app_llr[j][ldpc.row_edge[j][i]] / 2);
		}
		for (int i = 0; i < ldpc.row_w[j]; i++)
		{
			double temp = 1;
			for (int k = 0; k < ldpc.row_w[j]; k++)
			{
				if (k != i) temp *= ldpc.temp_llr[k];
			}
			if (temp == 1) ldpc.ext_llr[j][ldpc.row_edge[j][i]] = 40;
			else if (temp == 0) ldpc.ext_llr[j][ldpc.row_edge[j][i]] = 0;
			else if (temp == -1) ldpc.ext_llr[j][ldpc.row_edge[j][i]] = -40;
			else ldpc.ext_llr[j][ldpc.row_edge[j][i]] = log((1 + temp) / (1 - temp));
		}
	}
}

void LLR_V2C(LDPC &ldpc, LT &lt, int current_len)
{
	for (int i = 0; i < lt.in_len; i++)
	{
		double ldpc_temp = 0;
		for (int j = 0; j < ldpc.col_w[i]; j++)
		{
			ldpc_temp += ldpc.ext_llr[ldpc.col_edge[i][j]][i];
		}
		for (int j = 0; j < lt.i_deg[i]; j++)
		{
			if (lt.i_edge[i][j] < current_len)
			{
				double temp = ldpc_temp;
				for (int k = 0; k < lt.i_deg[i]; k++)
				{
					if ((k != j) && (lt.i_edge[i][k] < current_len)) temp += lt.llr_o2i[lt.i_edge[i][k]][lt.i_order[i][k]];
					if (lt.i_edge[i][k] >= current_len) break;
				}
				lt.llr_i2o[i][j] = temp;
			}
			else break;
		}
	}
	for (int i = 0; i < ldpc.h_col; i++)
	{
		double temp = 0;
		for (int j = 0; j < lt.i_deg[i]; j++)
		{
			if (lt.i_edge[i][j] < current_len) temp += lt.llr_o2i[lt.i_edge[i][j]][lt.i_order[i][j]];
		}
		for (int j = 0; j < ldpc.col_w[i]; j++)
		{
			ldpc.app_llr[ldpc.col_edge[i][j]][i] = temp;
			for (int k = 0; k < ldpc.col_w[i]; k++)
			{
				if (k != j) ldpc.app_llr[ldpc.col_edge[i][j]][i] += ldpc.ext_llr[ldpc.col_edge[i][k]][i];
			}
		}
	}
}