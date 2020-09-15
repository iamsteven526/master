#include <iostream>
#include "ldpc.h"
#include "lt.h"
#include "parameters.h"
using namespace std;

#pragma warning(disable:4996)

LDPC::LDPC(int column, int row)
{
	h_col = column;
	h_row = row;
	code_len = h_col;
	data_len = h_col - h_row;
	col_w = new unsigned short[h_col];				// weight of each column ( degrees of variable nodes )
	row_w = new unsigned short[h_row];				// weight of each row ( degrees of check nodes )
	col_edge = new unsigned short *[h_col];			// connected check nodes of each variable node 
	row_edge = new unsigned short *[h_row];			// connected variable nodes of each check node 
	G = new unsigned short *[h_col - h_row];
	unsigned short **H = new unsigned short *[h_row];
	unsigned short **H2 = new unsigned short *[h_row];
	for (int i = 0; i < h_row; i++)
	{
		H[i] = new unsigned short[h_col];
		H2[i] = new unsigned short[h_col];
		memset(H[i], 0, sizeof(unsigned short)*h_col);
		memset(H2[i], 0, sizeof(unsigned short)*h_col);
	}
	for (int i = 0; i < h_col - h_row; i++)
	{
		G[i] = new unsigned short[h_col];
		memset(G[i], 0, sizeof(unsigned short)*h_col);
	}
	ch_llr = new double[h_col];
	app_llr = new double *[h_row];
	ext_llr = new double *[h_row];
	for (int i = 0; i < h_row; i++)
	{
		app_llr[i] = new double[h_col];
		ext_llr[i] = new double[h_col];
	}
	//---------- read txt file ----------
	FILE *h_txt = fopen(INFILENAME, "r");
	for (int i = 0; i < h_col; i++)
	{
		fscanf(h_txt, " %hu", &col_w[i]);
	}
	for (int i = 0; i < h_row; i++)
	{
		fscanf(h_txt, " %hu", &row_w[i]);
	}
	for (int i = 0; i < h_col; i++)
	{
		col_edge[i] = new unsigned short[col_w[i]];
	}
	for (int i = 0; i < h_row; i++)
	{
		row_edge[i] = new unsigned short[row_w[i]];
	}
	for (int i = 0; i < h_col; i++)
	{
		int k = 0;
		for (int j = 0; j < col_w[i]; j++)
		{
			int temp;
			fscanf(h_txt, " %d", &temp);
			H[temp - 1][i] = 1;
			H2[temp - 1][i] = 1;
			col_edge[i][k++] = temp - 1;
		}
	}
	for (int i = 0; i < h_row; i++)
	{
		int k = 0;
		for (int j = 0; j < row_w[i]; j++)
		{
			int temp;
			fscanf(h_txt, " %d", &temp);
			row_edge[i][k++] = temp - 1;
		}
	}
	fclose(h_txt);
	//---------- Gaussian-Jordan elimination ----------
	for (int row = 0; row < h_row; row++)
	{
		if (H2[row][h_col - h_row + row] == 0)
		{
			for (int i = row + 1; i < h_row; i++)
			{
				if (H2[i][h_col - h_row + row] != 0)
				{
					for (int j = 0; j < h_col; j++)
					{
						unsigned short temp = H2[row][j];
						H2[row][j] = H2[i][j];
						H2[i][j] = temp;
					}
					break;
				}
			}
		}
		if (H2[row][h_col - h_row + row] == 0) // there's no rows have nonzero entry at the same column position
		{
			for (int i = 0; i < h_col; i++)
			{
				if (H2[row][i] != 0)
				{
					for (int j = 0; j < h_row; j++)
					{
						//---------- for row_edge ----------
						if (H[j][h_col - h_row + row] == 1)
						{
							for (int k = 0; k < row_w[j]; k++)
							{
								if (row_edge[j][k] == h_col - h_row + row)
								{
									row_edge[j][k] = i;
									break;
								}
							}
						}
						if (H[j][i] == 1)
						{
							for (int k = 0; k < row_w[j]; k++)
							{
								if (row_edge[j][k] == i)
								{
									row_edge[j][k] = h_col - h_row + row;
									break;
								}
							}
						}
						//---------- for column ----------
						unsigned short temp = H2[j][h_col - h_row + row];
						H2[j][h_col - h_row + row] = H2[j][i];
						H2[j][i] = temp;
						temp = H[j][h_col - h_row + row];
						H[j][h_col - h_row + row] = H[j][i];
						H[j][i] = temp;
					}
					//---------- for col_edge ----------
					Swap(&col_edge[h_col - h_row + row], &col_edge[i]);
					//---------- for col_w ----------
					unsigned short temp = col_w[h_col - h_row + row];
					col_w[h_col - h_row + row] = col_w[i];
					col_w[i] = temp;
					break;
				}
			}
		}
		if ((H2[row][h_col - h_row + row] == 1) && (row < h_row - 1))
		{
			for (int i = row + 1; i < h_row; i++)
			{
				if (H2[i][h_col - h_row + row] != 0)
				{
					for (int j = 0; j < h_col; j++)
					{
						H2[i][j] = H2[row][j] ^ H2[i][j];
					}
				}
			}
		}
	} // at the right-handed side, the bottom-left triangle entries are all zero now
	for (int row = h_row - 1; row > 0; row--)
	{
		if (H2[row][h_col - h_row + row] == 1)
		{
			for (int i = row - 1; i >= 0; i--)
			{
				if (H2[i][h_col - h_row + row] != 0)
				{
					for (int j = 0; j < h_col; j++)
					{
						H2[i][j] = H2[row][j] ^ H2[i][j];
					}
				}
			}
		}
	} // now H is the standard-form binary parity-check matrix
	//---------- generator matrix ----------
	for (int i = 0; i < h_col - h_row; i++)
	{
		G[i][i] = 1;
		for (int j = 0; j < h_row; j++)
		{
			G[i][h_col - h_row + j] = H2[j][i];
		}
	}
	//---------- check for the orthogonality ----------
	bool flag = true;
	int temp = 0;
	for (int i = 0; i < h_col - h_row; i++)
	{
		for (int j = 0; j < h_row; j++)
		{
			temp = 0;
			for (int k = 0; k < row_w[j]; k++)
			{
				temp = temp ^ G[i][row_edge[j][k]];
			}
			if (temp != 0) break;
		}
		if (temp != 0)
		{
			flag = false;
			break;
		}
	}
	if (!flag)
	{
		printf("Gaussian-Jordan elimination process is fail\n");
		system("pause");
	}
	for (int i = 0; i < h_row; i++)
	{
		delete[] H[i];
		delete[] H2[i];
	}
	delete[] H;
	delete[] H2;
}

void LDPC::Encoder(int *data, int *codeword)
{
	memcpy(codeword, data, sizeof(int)*data_len);
	for (int i = data_len; i < code_len; i++)
	{
		codeword[i] = 0;
		for (int j = 0; j < data_len; j++)
		{
			codeword[i] ^= data[j] * G[j][i];
		}
	}
}

void LDPC::Decoder(double variance, double *rx, double *decoded_llr, int it_max)
{
	int it = 0;
	//---------- initialization ----------
	for (int i = 0; i < h_col; i++)
	{
		ch_llr[i] = rx[i];
		for (int j = 0; j < h_row; j++)
		{
			app_llr[j][i] = ch_llr[i];
		}
	}
	while (true)
	{
		//---------- check messages ----------
		for (int j = 0; j < h_row; j++)
		{
			for (int i = 0; i < row_w[j]; i++)
			{
				if (app_llr[j][row_edge[j][i]] >= 30) { temp_llr[i] = 1;}
				else if (app_llr[j][row_edge[j][i]] <= -30) temp_llr[i] = -1;
				else temp_llr[i] = tanh(app_llr[j][row_edge[j][i]] / 2);
			}
			for (int i = 0; i < row_w[j]; i++)
			{
				double temp = 1;
				for (int k = 0; k < row_w[j]; k++)
				{
					if (k != i) temp *= temp_llr[k];
				}
				if (temp == 1) ext_llr[j][row_edge[j][i]] = 40;
				else if (temp == 0) ext_llr[j][row_edge[j][i]] = 0;
				else if (temp == -1) ext_llr[j][row_edge[j][i]] = -40;
				else ext_llr[j][row_edge[j][i]] = log((1 + temp) / (1 - temp));
			}
		}
		//---------- test ----------
		for (int i = 0; i < h_col; i++)
		{
			decoded_llr[i] = ch_llr[i];
			for (int j = 0; j < col_w[i]; j++)
			{
				decoded_llr[i] += ext_llr[col_edge[i][j]][i];
			}
		}
		bool flag = false;
		for (int i = 0; i < h_row; i++)
		{
			int temp = 0;
			for (int j = 0; j < row_w[i]; j++)
			{
				temp ^= HARD(decoded_llr[row_edge[i][j]]);
			}
			if (temp == 1)
			{
				flag = true;
				break;
			}
		}
		if (!flag || it == it_max) return; // early termination
		else
		{
			//---------- bit messages ----------
			for (int i = 0; i < h_col; i++)
			{
				for (int j = 0; j < col_w[i]; j++)
				{
					double temp = 0;
					for (int k = 0; k < col_w[i]; k++)
					{
						if (k != j) temp += ext_llr[col_edge[i][k]][i];
					}
					app_llr[col_edge[i][j]][i] = temp + ch_llr[i];
				}
			}
			it++;
		}
	}
}