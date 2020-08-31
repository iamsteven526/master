#include <iostream>
#include "ldpc.h"
#include "parameters.h"
#include <cmath>
#include <cstring>
using namespace std;
#pragma warning( disable : 4996 )

LDPC::LDPC(int column, int row)
{
	if (!CH_CODING_TYPE)
		return;
	h_col = column;											// number of coloums of parity-check matrix
	h_row = row;											// number of rows of parity-check matrix
	code_len = h_col;										// codeword length
	data_len = h_col - h_row;								// data length
	col_w = new unsigned short[h_col];						// weight of each column ( degree of variable node )
	row_w = new unsigned short[h_row];						// weight of each row ( degree of check node )
	col_edge = new unsigned short *[h_col];					// connected check nodes of each variable node 
	row_edge = new unsigned short *[h_row];					// connected variable nodes of each check node 
	G = new unsigned short *[h_col - h_row];				// generator matrix
	unsigned short **H = new unsigned short *[h_row];		// parity-check matrix
	unsigned short **H2 = new unsigned short *[h_row];		// parity-check matrix ( to performed Gaussain elimination )
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
	if (JCD)
	{
		p_v2c = new double **[h_col];
		p_c2v = new double **[h_row];
		for (int i = 0; i < h_col; i++)
		{
			p_v2c[i] = new double *[h_row];
		}
		for (int i = 0; i < h_row; i++)
		{
			p_c2v[i] = new double *[h_col];
		}
		for (int i = 0; i < h_col; i++)
		{
			for (int j = 0; j < h_row; j++)
			{
				p_v2c[i][j] = new double[NUM_LEVEL];
				p_c2v[j][i] = new double[NUM_LEVEL];
			}
		}
	}
	else
	{
		ch_llr = new double[h_col + JOINT_DEC*FFT_POINT];
		ch_value = new double[h_col + JOINT_DEC*FFT_POINT];
		app_llr = new double *[h_row];
		ext_llr = new double *[h_row];
		for (int i = 0; i < h_row; i++)
		{
			app_llr[i] = new double[h_col];
			ext_llr[i] = new double[h_col];
		}
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
		if (H2[row][h_col - h_row + row] == 0)
		{
			for (int i = 0; i < h_col; i++)
			{
				if (H2[row][i] != 0)
				{
					for (int j = 0; j < h_row; j++)
					{
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
						unsigned short temp = H2[j][h_col - h_row + row];
						H2[j][h_col - h_row + row] = H2[j][i];
						H2[j][i] = temp;
						temp = H[j][h_col - h_row + row];
						H[j][h_col - h_row + row] = H[j][i];
						H[j][i] = temp;
					}
					swap(&col_edge[h_col - h_row + row], &col_edge[i]);
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
	}
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
	}
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

	if (JOINT_DEC) DiffGraphConstruction();

}

void LDPC::Encoder(int *data, int *tx)
{
	memcpy(tx, data, sizeof(int)*data_len);
	for (int i = data_len; i < code_len; i++)
	{
		tx[i] = 0;
		for (int j = 0; j < data_len; j++)
		{
			tx[i] ^= data[j] * G[j][i];
		}
	}
}

void LDPC::SPA(double *rx, double *decoded_llr, int it_max)
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
		//---------- check-to-variable messages ----------
		for (int j = 0; j < h_row; j++)
		{
			for (int i = 0; i < row_w[j]; i++)
			{
				if (app_llr[j][row_edge[j][i]] >= 30) temp_value[i] = 1;
				else if (app_llr[j][row_edge[j][i]] <= -30) temp_value[i] = -1;
				else temp_value[i] = tanh(app_llr[j][row_edge[j][i]] / 2);
			}
			for (int i = 0; i < row_w[j]; i++)
			{
				double temp = 1;
				for (int k = 0; k < row_w[j]; k++)
				{
					if (k != i) temp *= temp_value[k];
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
				temp = temp ^ HARD(decoded_llr[row_edge[i][j]]);
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
			//---------- variable-to-check messages ----------
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

void LDPC::GSPA(double **app, int it_max)
{
	//---------- initialization ----------
	for (int i = 0; i < h_col; i++)
	{
		for (int j = 0; j < h_row; j++)
		{
			memcpy(p_v2c[i][j], app[i], sizeof(double)*NUM_LEVEL);
			memset(p_c2v[j][i], 0, sizeof(double)*NUM_LEVEL);
		}
	}
	//---------- iterative process ----------
	for (int it = 0; it < it_max; it++)
	{
		//---------- check-to-variable messages ----------
		for (int i = 0; i < h_row; i++)
		{
			for (int j = 0; j < row_w[i]; j++)
			{
				bool flag = true;
				for (int k = 0; k < row_w[i]; k++)
				{
					if (k != j)
					{
						if (flag)
						{
							memcpy(p_c2v[i][row_edge[i][j]], p_v2c[row_edge[i][k]][i], sizeof(double)*NUM_LEVEL);
							flag = false;
						}
						else
						{
							double temp[NUM_LEVEL] = { 0 };
							memcpy(temp, p_c2v[i][row_edge[i][j]], sizeof(double)*NUM_LEVEL);
							if (NUM_LEVEL == 4)
							{
								p_c2v[i][row_edge[i][j]][0] = temp[0] * p_v2c[row_edge[i][k]][i][0] + temp[1] * p_v2c[row_edge[i][k]][i][1] + temp[2] * p_v2c[row_edge[i][k]][i][2] + temp[3] * p_v2c[row_edge[i][k]][i][3];
								p_c2v[i][row_edge[i][j]][1] = temp[0] * p_v2c[row_edge[i][k]][i][1] + temp[1] * p_v2c[row_edge[i][k]][i][0] + temp[2] * p_v2c[row_edge[i][k]][i][3] + temp[3] * p_v2c[row_edge[i][k]][i][2];
								p_c2v[i][row_edge[i][j]][2] = temp[0] * p_v2c[row_edge[i][k]][i][2] + temp[1] * p_v2c[row_edge[i][k]][i][3] + temp[2] * p_v2c[row_edge[i][k]][i][0] + temp[3] * p_v2c[row_edge[i][k]][i][1];
								p_c2v[i][row_edge[i][j]][3] = temp[0] * p_v2c[row_edge[i][k]][i][3] + temp[1] * p_v2c[row_edge[i][k]][i][2] + temp[2] * p_v2c[row_edge[i][k]][i][1] + temp[3] * p_v2c[row_edge[i][k]][i][0];

							}
							else if (NUM_LEVEL == 8)
							{
								p_c2v[i][row_edge[i][j]][0] = temp[0] * p_v2c[row_edge[i][k]][i][0] + temp[1] * p_v2c[row_edge[i][k]][i][1] + temp[2] * p_v2c[row_edge[i][k]][i][2] + temp[3] * p_v2c[row_edge[i][k]][i][3] + temp[4] * p_v2c[row_edge[i][k]][i][4] + temp[5] * p_v2c[row_edge[i][k]][i][5] + temp[6] * p_v2c[row_edge[i][k]][i][6] + temp[7] * p_v2c[row_edge[i][k]][i][7];
								p_c2v[i][row_edge[i][j]][1] = temp[0] * p_v2c[row_edge[i][k]][i][1] + temp[1] * p_v2c[row_edge[i][k]][i][0] + temp[2] * p_v2c[row_edge[i][k]][i][3] + temp[3] * p_v2c[row_edge[i][k]][i][2] + temp[4] * p_v2c[row_edge[i][k]][i][5] + temp[5] * p_v2c[row_edge[i][k]][i][4] + temp[6] * p_v2c[row_edge[i][k]][i][7] + temp[7] * p_v2c[row_edge[i][k]][i][6];
								p_c2v[i][row_edge[i][j]][2] = temp[0] * p_v2c[row_edge[i][k]][i][2] + temp[1] * p_v2c[row_edge[i][k]][i][3] + temp[2] * p_v2c[row_edge[i][k]][i][0] + temp[3] * p_v2c[row_edge[i][k]][i][1] + temp[4] * p_v2c[row_edge[i][k]][i][6] + temp[5] * p_v2c[row_edge[i][k]][i][7] + temp[6] * p_v2c[row_edge[i][k]][i][4] + temp[7] * p_v2c[row_edge[i][k]][i][5];
								p_c2v[i][row_edge[i][j]][3] = temp[0] * p_v2c[row_edge[i][k]][i][3] + temp[1] * p_v2c[row_edge[i][k]][i][2] + temp[2] * p_v2c[row_edge[i][k]][i][1] + temp[3] * p_v2c[row_edge[i][k]][i][0] + temp[4] * p_v2c[row_edge[i][k]][i][7] + temp[5] * p_v2c[row_edge[i][k]][i][6] + temp[6] * p_v2c[row_edge[i][k]][i][5] + temp[7] * p_v2c[row_edge[i][k]][i][4];
								p_c2v[i][row_edge[i][j]][4] = temp[0] * p_v2c[row_edge[i][k]][i][4] + temp[1] * p_v2c[row_edge[i][k]][i][5] + temp[2] * p_v2c[row_edge[i][k]][i][6] + temp[3] * p_v2c[row_edge[i][k]][i][7] + temp[4] * p_v2c[row_edge[i][k]][i][0] + temp[5] * p_v2c[row_edge[i][k]][i][1] + temp[6] * p_v2c[row_edge[i][k]][i][2] + temp[7] * p_v2c[row_edge[i][k]][i][3];
								p_c2v[i][row_edge[i][j]][5] = temp[0] * p_v2c[row_edge[i][k]][i][5] + temp[1] * p_v2c[row_edge[i][k]][i][4] + temp[2] * p_v2c[row_edge[i][k]][i][7] + temp[3] * p_v2c[row_edge[i][k]][i][6] + temp[4] * p_v2c[row_edge[i][k]][i][1] + temp[5] * p_v2c[row_edge[i][k]][i][0] + temp[6] * p_v2c[row_edge[i][k]][i][3] + temp[7] * p_v2c[row_edge[i][k]][i][2];
								p_c2v[i][row_edge[i][j]][6] = temp[0] * p_v2c[row_edge[i][k]][i][6] + temp[1] * p_v2c[row_edge[i][k]][i][7] + temp[2] * p_v2c[row_edge[i][k]][i][4] + temp[3] * p_v2c[row_edge[i][k]][i][5] + temp[4] * p_v2c[row_edge[i][k]][i][2] + temp[5] * p_v2c[row_edge[i][k]][i][3] + temp[6] * p_v2c[row_edge[i][k]][i][0] + temp[7] * p_v2c[row_edge[i][k]][i][1];
								p_c2v[i][row_edge[i][j]][7] = temp[0] * p_v2c[row_edge[i][k]][i][7] + temp[1] * p_v2c[row_edge[i][k]][i][6] + temp[2] * p_v2c[row_edge[i][k]][i][5] + temp[3] * p_v2c[row_edge[i][k]][i][4] + temp[4] * p_v2c[row_edge[i][k]][i][3] + temp[5] * p_v2c[row_edge[i][k]][i][2] + temp[6] * p_v2c[row_edge[i][k]][i][1] + temp[7] * p_v2c[row_edge[i][k]][i][0];
							}
							else
							{
								printf("\nPARAMETER SETTING IS WRONG\n");
								system("pause");
							}
							for (int z = 0; z < NUM_LEVEL; z++)
							{
								if (p_c2v[i][row_edge[i][j]][z] < NUMERIC_LIMIT) p_c2v[i][j][z] = NUMERIC_LIMIT;
							}
						}
					}
				}
			}
		}
		//---------- variable-to-check messages ----------
		for (int i = 0; i < h_col; i++)
		{
			for (int j = 0; j < col_w[i]; j++)
			{
				memcpy(p_v2c[i][col_edge[i][j]], app[i], sizeof(double)*NUM_LEVEL);
				for (int k = 0; k < col_w[i]; k++)
				{
					if (k != j)
					{
						for (int z = 0; z < NUM_LEVEL; z++)
						{
							p_v2c[i][col_edge[i][j]][z] *= p_c2v[col_edge[i][k]][i][z];
						}
						double temp = 0;
						for (int z = 0; z < NUM_LEVEL; z++)
						{
							if (p_v2c[i][col_edge[i][j]][z] < NUMERIC_LIMIT) p_v2c[i][col_edge[i][j]][z] = NUMERIC_LIMIT;
							temp += p_v2c[i][col_edge[i][j]][z];
						}
						for (int z = 0; z < NUM_LEVEL; z++)
						{
							p_v2c[i][col_edge[i][j]][z] /= temp;
							if (p_v2c[i][col_edge[i][j]][z] < NUMERIC_LIMIT) p_v2c[i][col_edge[i][j]][z] = NUMERIC_LIMIT;
						}
					}
				}
			}
		}
	}
	//---------- finalization ---------- 
	for (int i = 0; i < h_col; i++)
	{
		for (int j = 0; j < col_w[i]; j++)
		{
			for (int z = 0; z < NUM_LEVEL; z++)
			{
				app[i][z] *= p_c2v[col_edge[i][j]][i][z];
			}
			double temp = 0;
			for (int z = 0; z < NUM_LEVEL; z++)
			{
				if (app[i][z] < NUMERIC_LIMIT) app[i][z] = NUMERIC_LIMIT;
				temp += app[i][z];
			}
			for (int z = 0; z < NUM_LEVEL; z++)
			{
				app[i][z] /= temp;
				if (app[i][z] < NUMERIC_LIMIT) app[i][z] = NUMERIC_LIMIT;
			}
		}
	}
}

void LDPC::DiffGraphConstruction()
{
	out_len = code_len + FFT_POINT;
	o_deg = new int[out_len];
	i_deg = new int[code_len];
	o_edge = new int*[out_len];
	i_edge = new int*[code_len];
	o_link = new int*[out_len];
	i_link = new int*[code_len];
	for (int i = 0; i < code_len; i++)
	{
		i_deg[i] = 2;
		i_edge[i] = new int[i_deg[i]];
		i_link[i] = new int[i_deg[i]];
	}
	for (int i = 0; i < FFT_POINT; i++)
	{
		o_deg[i] = 1;
		o_edge[i] = new int[o_deg[i]];
		o_link[i] = new int[o_deg[i]];
		o_deg[out_len - i - 1] = 1;
		o_edge[out_len - i - 1] = new int[o_deg[out_len - i]];
		o_link[out_len - i - 1] = new int[o_deg[out_len - i]];
	}
	for (int i = FFT_POINT; i < out_len - FFT_POINT; i++)
	{
		o_deg[i] = 2;
		o_edge[i] = new int[o_deg[i]];
		o_link[i] = new int[o_deg[i]];
	}
	int *count = new int[out_len];
	memset(count, 0, sizeof(int)*out_len);
	for (int i = 0; i < FFT_POINT; i++)
	{
		for (int j = 0; j < FFT_SEGMENT; j++)
		{
			i_edge[j*FFT_POINT + i][0] = j*FFT_POINT + i;
			i_link[j*FFT_POINT + i][0] = count[j*FFT_POINT + i];
			o_edge[j*FFT_POINT + i][count[j*FFT_POINT + i]] = j*FFT_POINT + i;
			o_link[j*FFT_POINT + i][count[j*FFT_POINT + i]] = 0;
			count[j*FFT_POINT + i]++;
			i_edge[j*FFT_POINT + i][1] = (j + 1)*FFT_POINT + i;
			i_link[j*FFT_POINT + i][1] = count[(j + 1)*FFT_POINT + i];
			o_edge[(j + 1)*FFT_POINT + i][count[(j + 1)*FFT_POINT + i]] = j*FFT_POINT + i;
			o_link[(j + 1)*FFT_POINT + i][count[(j + 1)*FFT_POINT + i]] = 1;
			count[(j + 1)*FFT_POINT + i]++;
		}
	}
	delete[] count;
	if (JCD)
	{
		p_i2o = new double **[h_col];
		for (int i = 0; i < h_col; i++)
		{
			p_i2o[i] = new double *[2];
			for (int j = 0; j < 2; j++)
			{
				p_i2o[i][j] = new double[NUM_LEVEL];
			}
		}
		p_o2i = new double **[out_len];
		for (int i = 0; i < FFT_POINT; i++)
		{
			p_o2i[i] = new double *[1];
			p_o2i[out_len - i - 1] = new double *[1];
			for (int j = 0; j < 1; j++)
			{
				p_o2i[i][j] = new double[NUM_LEVEL];
				p_o2i[out_len - i - 1][j] = new double[NUM_LEVEL];
			}
		}
		for (int i = FFT_POINT; i < out_len - FFT_POINT; i++)
		{
			p_o2i[i] = new double *[2];
			for (int j = 0; j < 2; j++)
			{
				p_o2i[i][j] = new double[NUM_LEVEL];
			}
		}
		input_p = new double*[h_col];
		output_p = new double*[h_col];
		for (int i = 0; i < h_col; i++)
		{
			input_p[i] = new double[NUM_LEVEL];
			output_p[i] = new double[NUM_LEVEL];
		}
	}
	else
	{
		llr_o2i = new double*[out_len];
		llr_i2o = new double*[code_len];
		for (int i = 0; i < code_len; i++)
		{
			llr_i2o[i] = new double[i_deg[i]];
		}
		for (int i = 0; i < FFT_POINT; i++)
		{
			llr_o2i[i] = new double[o_deg[i]];
			llr_o2i[out_len - i - 1] = new double[o_deg[out_len - i]];
		}
		for (int i = FFT_POINT; i < out_len - FFT_POINT; i++)
		{
			llr_o2i[i] = new double[o_deg[i]];
		}
		input_llr = new double[code_len];
		output_llr = new double[code_len];
	}
}

void LDPC::JointO2I()
{
	for (int i = 0; i < out_len; i++)
	{
		if (o_deg[i] == 1)
		{
			llr_o2i[i][0] = ch_llr[i];
		}
		else
		{
			for (int j = 0; j < o_deg[i]; j++)
			{
				llr_o2i[i][j] = ch_llr[i];
				for (int k = 0; k < o_deg[i]; k++)
				{
					if (k != j) llr_o2i[i][j] += llr_i2o[o_edge[i][k]][o_link[i][k]];
				}
			}
		}
	}
}

void LDPC::JointI2O(double *input, double *output)
{
	for (int i = 0; i < code_len; i++)
	{
		for (int j = 0; j < i_deg[i]; j++)
		{
			if (llr_o2i[i_edge[i][j]][i_link[i][j]] >= 30) temp_value[j] = 1;
			else if (llr_o2i[i_edge[i][j]][i_link[i][j]] <= -30) temp_value[j] = -1;
			else temp_value[j] = tanh(llr_o2i[i_edge[i][j]][i_link[i][j]] / 2);
		}
		for (int j = 0; j < i_deg[i]; j++)
		{
			double temp = tanh(input[i] / 2);
			for (int k = 0; k < i_deg[i]; k++)
			{
				if (k != j) temp *= temp_value[k];
			}
			if (temp == 1) llr_i2o[i][j] = 40;
			else if (temp == 0) llr_i2o[i][j] = 0;
			else if (temp == -1) llr_i2o[i][j] = -40;
			else llr_i2o[i][j] = log((1 + temp) / (1 - temp));
		}
		output[i] = 1;
		for (int j = 0; j < i_deg[i]; j++)
		{
			output[i] *= temp_value[j];
		}
		if (output[i] == 1) output[i] = 40;
		else if (output[i] == 0) output[i] = 0;
		else if (output[i] == -1) output[i] = -40;
		else output[i] = log((1 + output[i]) / (1 - output[i]));
	}
}

void LDPC::JointSPA(double *appLlr, double *refLlr, double *decodedLlr, int outer_it_max, int inter_it_max)
{
	//---------- initialization ----------
	for (int i = 0; i < h_row; i++)
	{
		memset(app_llr[i], 0, h_col * sizeof(double));
	}
	for (int i = 0; i < (out_len - code_len); i++)
	{
		ch_llr[i] = refLlr[i];
		if (ch_llr[i] >= 30) ch_value[i] = 1;
		else if (ch_llr[i] <= -30) ch_value[i] = -1;
		else ch_value[i] = tanh(ch_llr[i] / 2);
	}
	for (int i = (out_len - code_len); i < out_len; i++)
	{
		ch_llr[i] = appLlr[i - (out_len - code_len)];
		if (ch_llr[i] >= 30) ch_value[i] = 1;
		else if (ch_llr[i] <= -30) ch_value[i] = -1;
		else ch_value[i] = tanh(ch_llr[i] / 2);
	}
	for (int i = 0; i < out_len; i++)
	{
		memset(llr_o2i[i], 0, sizeof(double)*o_deg[i]);
	}
	for (int i = 0; i < code_len; i++)
	{
		memset(llr_i2o[i], 0, sizeof(double)*i_deg[i]);
	}
	memset(input_llr, 0, sizeof(double)*code_len);
	memset(output_llr, 0, sizeof(double)*code_len);
	//---------- belief propagation ----------
	for (int it = 1; it <= outer_it_max; it++)
	{
		//---------- differential component ---------- 
		for (int inter_it = 1; inter_it <= inter_it_max; inter_it++)
		{
			JointO2I();
			JointI2O(input_llr, output_llr);
		}
		//---------- LDPC component ---------- 
		for (int inter_it = 1; inter_it <= inter_it_max; inter_it++)
		{
			for (int j = 0; j < h_row; j++)
			{
				for (int i = 0; i < row_w[j]; i++)
				{
					if (app_llr[j][row_edge[j][i]] >= 30) temp_value[i] = 1;
					else if (app_llr[j][row_edge[j][i]] <= -30) temp_value[i] = -1;
					else temp_value[i] = tanh(app_llr[j][row_edge[j][i]] / 2);
				}
				for (int i = 0; i < row_w[j]; i++)
				{
					double temp = 1;
					for (int k = 0; k < row_w[j]; k++)
					{
						if (k != i) temp *= temp_value[k];
					}
					if (temp == 1) ext_llr[j][row_edge[j][i]] = 40;
					else if (temp == 0) ext_llr[j][row_edge[j][i]] = 0;
					else if (temp == -1) ext_llr[j][row_edge[j][i]] = -40;
					else ext_llr[j][row_edge[j][i]] = log((1 + temp) / (1 - temp));
				}
			}
			for (int i = 0; i < code_len; i++)
			{
				decodedLlr[i] = output_llr[i];
				for (int j = 0; j < col_w[i]; j++)
				{
					decodedLlr[i] += ext_llr[col_edge[i][j]][i];
				}
			}
			if (it > 5)
			{
				bool flag = false;
				for (int i = 0; i < h_row; i++)
				{
					int temp = 0;
					for (int j = 0; j < row_w[i]; j++)
					{
						temp ^= HARD(decodedLlr[row_edge[i][j]]);
					}
					if (temp == 1)
					{
						flag = true;
						break;
					}
				}
				if (!flag) break; // early termination
			}
			for (int i = 0; i < h_col; i++)
			{
				for (int j = 0; j < col_w[i]; j++)
				{
					app_llr[col_edge[i][j]][i] = output_llr[i];
					for (int k = 0; k < col_w[i]; k++)
					{
						if (k != j) app_llr[col_edge[i][j]][i] += ext_llr[col_edge[i][k]][i];
					}
				}
			}
			for (int i = 0; i < h_col; i++)
			{
				input_llr[i] = 0;
				for (int j = 0; j < col_w[i]; j++)
				{
					input_llr[i] += ext_llr[col_edge[i][j]][i];
				}
			}
		}
	}
}

void LDPC::JointGSPA(double **app, int outer_it_max, int inter_it_max)
{
	//---------- initialization ----------
	for (int i = 0; i < h_col; i++)
	{
		for (int j = 0; j < h_row; j++)
		{
			memset(p_v2c[i][j], 0, sizeof(double)*NUM_LEVEL);
			memset(p_c2v[j][i], 0, sizeof(double)*NUM_LEVEL);
		}
		memset(input_p[i], 0, sizeof(double)*NUM_LEVEL);
		memset(output_p[i], 0, sizeof(double)*NUM_LEVEL);
		for (int j = 0; j < 2; j++)
		{
			memset(p_i2o[i][j], 0, sizeof(double)*NUM_LEVEL);
		}
	}
	for (int i = 0; i < FFT_POINT; i++)
	{
		memset(p_o2i[i][0], 0, sizeof(double)*NUM_LEVEL);
		memset(p_o2i[out_len - i - 1][0], 0, sizeof(double)*NUM_LEVEL);
	}
	for (int i = FFT_POINT; i < out_len - FFT_POINT; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			memset(p_o2i[i][j], 0, sizeof(double)*NUM_LEVEL);
		}
	}
	//---------- iterative process ----------
	for (int it = 1; it <= outer_it_max; it++)
	{
		//---------- differential component ---------- 
		for (int inter_it = 1; inter_it <= inter_it_max; inter_it++)
		{
			//---------- variable-to-check messages ----------
			for (int i = 0; i < out_len; i++)
			{
				for (int j = 0; j < o_deg[i]; j++)
				{
					memcpy(p_o2i[i][j], app[i], sizeof(double)*NUM_LEVEL);
					for (int k = 0; k < o_deg[i]; k++)
					{
						if (k != j)
						{
							for (int z = 0; z < NUM_LEVEL; z++)
							{
								p_o2i[i][j][z] *= p_i2o[o_edge[i][k]][o_link[i][k]][z];
							}
							double temp = 0;
							for (int z = 0; z < NUM_LEVEL; z++)
							{
								if (p_o2i[i][j][z] < NUMERIC_LIMIT) p_o2i[i][j][z] = NUMERIC_LIMIT;
								temp += p_o2i[i][j][z];
							}
							for (int z = 0; z < NUM_LEVEL; z++)
							{
								p_o2i[i][j][z] /= temp;
								if (p_o2i[i][j][z] < NUMERIC_LIMIT) p_o2i[i][j][z] = NUMERIC_LIMIT;
							}
						}
					}
				}
			}
			//---------- check-to-variable messages ----------
			for (int i = 0; i < h_col; i++)
			{
				for (int j = 0; j < i_deg[i]; j++)
				{
					for (int k = 0; k < i_deg[i]; k++)
					{
						if (k != j)
						{
							if (NUM_LEVEL == 4)
							{
								p_i2o[i][j][0] = input_p[i][0] * p_o2i[i_edge[i][k]][i_link[i][k]][0] + input_p[i][1] * p_o2i[i_edge[i][k]][i_link[i][k]][1] + input_p[i][2] * p_o2i[i_edge[i][k]][i_link[i][k]][2] + input_p[i][3] * p_o2i[i_edge[i][k]][i_link[i][k]][3];
								p_i2o[i][j][1] = input_p[i][0] * p_o2i[i_edge[i][k]][i_link[i][k]][1] + input_p[i][1] * p_o2i[i_edge[i][k]][i_link[i][k]][0] + input_p[i][2] * p_o2i[i_edge[i][k]][i_link[i][k]][3] + input_p[i][3] * p_o2i[i_edge[i][k]][i_link[i][k]][2];
								p_i2o[i][j][2] = input_p[i][0] * p_o2i[i_edge[i][k]][i_link[i][k]][2] + input_p[i][1] * p_o2i[i_edge[i][k]][i_link[i][k]][3] + input_p[i][2] * p_o2i[i_edge[i][k]][i_link[i][k]][0] + input_p[i][3] * p_o2i[i_edge[i][k]][i_link[i][k]][1];
								p_i2o[i][j][3] = input_p[i][0] * p_o2i[i_edge[i][k]][i_link[i][k]][3] + input_p[i][1] * p_o2i[i_edge[i][k]][i_link[i][k]][2] + input_p[i][2] * p_o2i[i_edge[i][k]][i_link[i][k]][1] + input_p[i][3] * p_o2i[i_edge[i][k]][i_link[i][k]][0];
							}
							else if (NUM_LEVEL == 8)
							{
								p_i2o[i][j][0] = input_p[i][0] * p_o2i[i_edge[i][k]][i_link[i][k]][0] + input_p[i][1] * p_o2i[i_edge[i][k]][i_link[i][k]][1] + input_p[i][2] * p_o2i[i_edge[i][k]][i_link[i][k]][2] + input_p[i][3] * p_o2i[i_edge[i][k]][i_link[i][k]][3] + input_p[i][4] * p_o2i[i_edge[i][k]][i_link[i][k]][4] + input_p[i][5] * p_o2i[i_edge[i][k]][i_link[i][k]][5] + input_p[i][6] * p_o2i[i_edge[i][k]][i_link[i][k]][6] + input_p[i][7] * p_o2i[i_edge[i][k]][i_link[i][k]][7];
								p_i2o[i][j][1] = input_p[i][0] * p_o2i[i_edge[i][k]][i_link[i][k]][1] + input_p[i][1] * p_o2i[i_edge[i][k]][i_link[i][k]][0] + input_p[i][2] * p_o2i[i_edge[i][k]][i_link[i][k]][3] + input_p[i][3] * p_o2i[i_edge[i][k]][i_link[i][k]][2] + input_p[i][4] * p_o2i[i_edge[i][k]][i_link[i][k]][5] + input_p[i][5] * p_o2i[i_edge[i][k]][i_link[i][k]][4] + input_p[i][6] * p_o2i[i_edge[i][k]][i_link[i][k]][7] + input_p[i][7] * p_o2i[i_edge[i][k]][i_link[i][k]][6];
								p_i2o[i][j][2] = input_p[i][0] * p_o2i[i_edge[i][k]][i_link[i][k]][2] + input_p[i][1] * p_o2i[i_edge[i][k]][i_link[i][k]][3] + input_p[i][2] * p_o2i[i_edge[i][k]][i_link[i][k]][0] + input_p[i][3] * p_o2i[i_edge[i][k]][i_link[i][k]][1] + input_p[i][4] * p_o2i[i_edge[i][k]][i_link[i][k]][6] + input_p[i][5] * p_o2i[i_edge[i][k]][i_link[i][k]][7] + input_p[i][6] * p_o2i[i_edge[i][k]][i_link[i][k]][4] + input_p[i][7] * p_o2i[i_edge[i][k]][i_link[i][k]][5];
								p_i2o[i][j][3] = input_p[i][0] * p_o2i[i_edge[i][k]][i_link[i][k]][3] + input_p[i][1] * p_o2i[i_edge[i][k]][i_link[i][k]][2] + input_p[i][2] * p_o2i[i_edge[i][k]][i_link[i][k]][1] + input_p[i][3] * p_o2i[i_edge[i][k]][i_link[i][k]][0] + input_p[i][4] * p_o2i[i_edge[i][k]][i_link[i][k]][7] + input_p[i][5] * p_o2i[i_edge[i][k]][i_link[i][k]][6] + input_p[i][6] * p_o2i[i_edge[i][k]][i_link[i][k]][5] + input_p[i][7] * p_o2i[i_edge[i][k]][i_link[i][k]][4];
								p_i2o[i][j][4] = input_p[i][0] * p_o2i[i_edge[i][k]][i_link[i][k]][4] + input_p[i][1] * p_o2i[i_edge[i][k]][i_link[i][k]][5] + input_p[i][2] * p_o2i[i_edge[i][k]][i_link[i][k]][6] + input_p[i][3] * p_o2i[i_edge[i][k]][i_link[i][k]][7] + input_p[i][4] * p_o2i[i_edge[i][k]][i_link[i][k]][0] + input_p[i][5] * p_o2i[i_edge[i][k]][i_link[i][k]][1] + input_p[i][6] * p_o2i[i_edge[i][k]][i_link[i][k]][2] + input_p[i][7] * p_o2i[i_edge[i][k]][i_link[i][k]][3];
								p_i2o[i][j][5] = input_p[i][0] * p_o2i[i_edge[i][k]][i_link[i][k]][5] + input_p[i][1] * p_o2i[i_edge[i][k]][i_link[i][k]][4] + input_p[i][2] * p_o2i[i_edge[i][k]][i_link[i][k]][7] + input_p[i][3] * p_o2i[i_edge[i][k]][i_link[i][k]][6] + input_p[i][4] * p_o2i[i_edge[i][k]][i_link[i][k]][1] + input_p[i][5] * p_o2i[i_edge[i][k]][i_link[i][k]][0] + input_p[i][6] * p_o2i[i_edge[i][k]][i_link[i][k]][3] + input_p[i][7] * p_o2i[i_edge[i][k]][i_link[i][k]][2];
								p_i2o[i][j][6] = input_p[i][0] * p_o2i[i_edge[i][k]][i_link[i][k]][6] + input_p[i][1] * p_o2i[i_edge[i][k]][i_link[i][k]][7] + input_p[i][2] * p_o2i[i_edge[i][k]][i_link[i][k]][4] + input_p[i][3] * p_o2i[i_edge[i][k]][i_link[i][k]][5] + input_p[i][4] * p_o2i[i_edge[i][k]][i_link[i][k]][2] + input_p[i][5] * p_o2i[i_edge[i][k]][i_link[i][k]][3] + input_p[i][6] * p_o2i[i_edge[i][k]][i_link[i][k]][0] + input_p[i][7] * p_o2i[i_edge[i][k]][i_link[i][k]][1];
								p_i2o[i][j][7] = input_p[i][0] * p_o2i[i_edge[i][k]][i_link[i][k]][7] + input_p[i][1] * p_o2i[i_edge[i][k]][i_link[i][k]][6] + input_p[i][2] * p_o2i[i_edge[i][k]][i_link[i][k]][5] + input_p[i][3] * p_o2i[i_edge[i][k]][i_link[i][k]][4] + input_p[i][4] * p_o2i[i_edge[i][k]][i_link[i][k]][3] + input_p[i][5] * p_o2i[i_edge[i][k]][i_link[i][k]][2] + input_p[i][6] * p_o2i[i_edge[i][k]][i_link[i][k]][1] + input_p[i][7] * p_o2i[i_edge[i][k]][i_link[i][k]][0];
							}
							else
							{
								printf("\nPARAMETER SETTING IS WRONG\n");
								system("pause");
							}
							for (int z = 0; z < NUM_LEVEL; z++)
							{
								if (p_i2o[i][j][z] < NUMERIC_LIMIT) p_i2o[i][j][z] = NUMERIC_LIMIT;
							}
						}
					}
				}
				if (NUM_LEVEL == 4)
				{
					output_p[i][0] = p_o2i[i_edge[i][0]][i_link[i][0]][0] * p_o2i[i_edge[i][1]][i_link[i][1]][0] + p_o2i[i_edge[i][0]][i_link[i][0]][1] * p_o2i[i_edge[i][1]][i_link[i][1]][1] + p_o2i[i_edge[i][0]][i_link[i][0]][2] * p_o2i[i_edge[i][1]][i_link[i][1]][2] + p_o2i[i_edge[i][0]][i_link[i][0]][3] * p_o2i[i_edge[i][1]][i_link[i][1]][3];
					output_p[i][1] = p_o2i[i_edge[i][0]][i_link[i][0]][0] * p_o2i[i_edge[i][1]][i_link[i][1]][1] + p_o2i[i_edge[i][0]][i_link[i][0]][1] * p_o2i[i_edge[i][1]][i_link[i][1]][0] + p_o2i[i_edge[i][0]][i_link[i][0]][2] * p_o2i[i_edge[i][1]][i_link[i][1]][3] + p_o2i[i_edge[i][0]][i_link[i][0]][3] * p_o2i[i_edge[i][1]][i_link[i][1]][2];
					output_p[i][2] = p_o2i[i_edge[i][0]][i_link[i][0]][0] * p_o2i[i_edge[i][1]][i_link[i][1]][2] + p_o2i[i_edge[i][0]][i_link[i][0]][1] * p_o2i[i_edge[i][1]][i_link[i][1]][3] + p_o2i[i_edge[i][0]][i_link[i][0]][2] * p_o2i[i_edge[i][1]][i_link[i][1]][0] + p_o2i[i_edge[i][0]][i_link[i][0]][3] * p_o2i[i_edge[i][1]][i_link[i][1]][1];
					output_p[i][3] = p_o2i[i_edge[i][0]][i_link[i][0]][0] * p_o2i[i_edge[i][1]][i_link[i][1]][3] + p_o2i[i_edge[i][0]][i_link[i][0]][1] * p_o2i[i_edge[i][1]][i_link[i][1]][2] + p_o2i[i_edge[i][0]][i_link[i][0]][2] * p_o2i[i_edge[i][1]][i_link[i][1]][1] + p_o2i[i_edge[i][0]][i_link[i][0]][3] * p_o2i[i_edge[i][1]][i_link[i][1]][0];
				}
				else if (NUM_LEVEL == 8)
				{
					output_p[i][0] = p_o2i[i_edge[i][0]][i_link[i][0]][0] * p_o2i[i_edge[i][1]][i_link[i][1]][0] + p_o2i[i_edge[i][0]][i_link[i][0]][1] * p_o2i[i_edge[i][1]][i_link[i][1]][1] + p_o2i[i_edge[i][0]][i_link[i][0]][2] * p_o2i[i_edge[i][1]][i_link[i][1]][2] + p_o2i[i_edge[i][0]][i_link[i][0]][3] * p_o2i[i_edge[i][1]][i_link[i][1]][3] + p_o2i[i_edge[i][0]][i_link[i][0]][4] * p_o2i[i_edge[i][1]][i_link[i][1]][4] + p_o2i[i_edge[i][0]][i_link[i][0]][5] * p_o2i[i_edge[i][1]][i_link[i][1]][5] + p_o2i[i_edge[i][0]][i_link[i][0]][6] * p_o2i[i_edge[i][1]][i_link[i][1]][6] + p_o2i[i_edge[i][0]][i_link[i][0]][7] * p_o2i[i_edge[i][1]][i_link[i][1]][7];
					output_p[i][1] = p_o2i[i_edge[i][0]][i_link[i][0]][0] * p_o2i[i_edge[i][1]][i_link[i][1]][1] + p_o2i[i_edge[i][0]][i_link[i][0]][1] * p_o2i[i_edge[i][1]][i_link[i][1]][0] + p_o2i[i_edge[i][0]][i_link[i][0]][2] * p_o2i[i_edge[i][1]][i_link[i][1]][3] + p_o2i[i_edge[i][0]][i_link[i][0]][3] * p_o2i[i_edge[i][1]][i_link[i][1]][2] + p_o2i[i_edge[i][0]][i_link[i][0]][4] * p_o2i[i_edge[i][1]][i_link[i][1]][5] + p_o2i[i_edge[i][0]][i_link[i][0]][5] * p_o2i[i_edge[i][1]][i_link[i][1]][4] + p_o2i[i_edge[i][0]][i_link[i][0]][6] * p_o2i[i_edge[i][1]][i_link[i][1]][7] + p_o2i[i_edge[i][0]][i_link[i][0]][7] * p_o2i[i_edge[i][1]][i_link[i][1]][6];
					output_p[i][2] = p_o2i[i_edge[i][0]][i_link[i][0]][0] * p_o2i[i_edge[i][1]][i_link[i][1]][2] + p_o2i[i_edge[i][0]][i_link[i][0]][1] * p_o2i[i_edge[i][1]][i_link[i][1]][3] + p_o2i[i_edge[i][0]][i_link[i][0]][2] * p_o2i[i_edge[i][1]][i_link[i][1]][0] + p_o2i[i_edge[i][0]][i_link[i][0]][3] * p_o2i[i_edge[i][1]][i_link[i][1]][1] + p_o2i[i_edge[i][0]][i_link[i][0]][4] * p_o2i[i_edge[i][1]][i_link[i][1]][6] + p_o2i[i_edge[i][0]][i_link[i][0]][5] * p_o2i[i_edge[i][1]][i_link[i][1]][7] + p_o2i[i_edge[i][0]][i_link[i][0]][6] * p_o2i[i_edge[i][1]][i_link[i][1]][4] + p_o2i[i_edge[i][0]][i_link[i][0]][7] * p_o2i[i_edge[i][1]][i_link[i][1]][5];
					output_p[i][3] = p_o2i[i_edge[i][0]][i_link[i][0]][0] * p_o2i[i_edge[i][1]][i_link[i][1]][3] + p_o2i[i_edge[i][0]][i_link[i][0]][1] * p_o2i[i_edge[i][1]][i_link[i][1]][2] + p_o2i[i_edge[i][0]][i_link[i][0]][2] * p_o2i[i_edge[i][1]][i_link[i][1]][1] + p_o2i[i_edge[i][0]][i_link[i][0]][3] * p_o2i[i_edge[i][1]][i_link[i][1]][0] + p_o2i[i_edge[i][0]][i_link[i][0]][4] * p_o2i[i_edge[i][1]][i_link[i][1]][7] + p_o2i[i_edge[i][0]][i_link[i][0]][5] * p_o2i[i_edge[i][1]][i_link[i][1]][6] + p_o2i[i_edge[i][0]][i_link[i][0]][6] * p_o2i[i_edge[i][1]][i_link[i][1]][5] + p_o2i[i_edge[i][0]][i_link[i][0]][7] * p_o2i[i_edge[i][1]][i_link[i][1]][4];
					output_p[i][4] = p_o2i[i_edge[i][0]][i_link[i][0]][0] * p_o2i[i_edge[i][1]][i_link[i][1]][4] + p_o2i[i_edge[i][0]][i_link[i][0]][1] * p_o2i[i_edge[i][1]][i_link[i][1]][5] + p_o2i[i_edge[i][0]][i_link[i][0]][2] * p_o2i[i_edge[i][1]][i_link[i][1]][6] + p_o2i[i_edge[i][0]][i_link[i][0]][3] * p_o2i[i_edge[i][1]][i_link[i][1]][7] + p_o2i[i_edge[i][0]][i_link[i][0]][4] * p_o2i[i_edge[i][1]][i_link[i][1]][0] + p_o2i[i_edge[i][0]][i_link[i][0]][5] * p_o2i[i_edge[i][1]][i_link[i][1]][1] + p_o2i[i_edge[i][0]][i_link[i][0]][6] * p_o2i[i_edge[i][1]][i_link[i][1]][2] + p_o2i[i_edge[i][0]][i_link[i][0]][7] * p_o2i[i_edge[i][1]][i_link[i][1]][3];
					output_p[i][5] = p_o2i[i_edge[i][0]][i_link[i][0]][0] * p_o2i[i_edge[i][1]][i_link[i][1]][5] + p_o2i[i_edge[i][0]][i_link[i][0]][1] * p_o2i[i_edge[i][1]][i_link[i][1]][4] + p_o2i[i_edge[i][0]][i_link[i][0]][2] * p_o2i[i_edge[i][1]][i_link[i][1]][7] + p_o2i[i_edge[i][0]][i_link[i][0]][3] * p_o2i[i_edge[i][1]][i_link[i][1]][6] + p_o2i[i_edge[i][0]][i_link[i][0]][4] * p_o2i[i_edge[i][1]][i_link[i][1]][1] + p_o2i[i_edge[i][0]][i_link[i][0]][5] * p_o2i[i_edge[i][1]][i_link[i][1]][0] + p_o2i[i_edge[i][0]][i_link[i][0]][6] * p_o2i[i_edge[i][1]][i_link[i][1]][3] + p_o2i[i_edge[i][0]][i_link[i][0]][7] * p_o2i[i_edge[i][1]][i_link[i][1]][2];
					output_p[i][6] = p_o2i[i_edge[i][0]][i_link[i][0]][0] * p_o2i[i_edge[i][1]][i_link[i][1]][6] + p_o2i[i_edge[i][0]][i_link[i][0]][1] * p_o2i[i_edge[i][1]][i_link[i][1]][7] + p_o2i[i_edge[i][0]][i_link[i][0]][2] * p_o2i[i_edge[i][1]][i_link[i][1]][4] + p_o2i[i_edge[i][0]][i_link[i][0]][3] * p_o2i[i_edge[i][1]][i_link[i][1]][5] + p_o2i[i_edge[i][0]][i_link[i][0]][4] * p_o2i[i_edge[i][1]][i_link[i][1]][2] + p_o2i[i_edge[i][0]][i_link[i][0]][5] * p_o2i[i_edge[i][1]][i_link[i][1]][3] + p_o2i[i_edge[i][0]][i_link[i][0]][6] * p_o2i[i_edge[i][1]][i_link[i][1]][0] + p_o2i[i_edge[i][0]][i_link[i][0]][7] * p_o2i[i_edge[i][1]][i_link[i][1]][1];
					output_p[i][7] = p_o2i[i_edge[i][0]][i_link[i][0]][0] * p_o2i[i_edge[i][1]][i_link[i][1]][7] + p_o2i[i_edge[i][0]][i_link[i][0]][1] * p_o2i[i_edge[i][1]][i_link[i][1]][6] + p_o2i[i_edge[i][0]][i_link[i][0]][2] * p_o2i[i_edge[i][1]][i_link[i][1]][5] + p_o2i[i_edge[i][0]][i_link[i][0]][3] * p_o2i[i_edge[i][1]][i_link[i][1]][4] + p_o2i[i_edge[i][0]][i_link[i][0]][4] * p_o2i[i_edge[i][1]][i_link[i][1]][3] + p_o2i[i_edge[i][0]][i_link[i][0]][5] * p_o2i[i_edge[i][1]][i_link[i][1]][2] + p_o2i[i_edge[i][0]][i_link[i][0]][6] * p_o2i[i_edge[i][1]][i_link[i][1]][1] + p_o2i[i_edge[i][0]][i_link[i][0]][7] * p_o2i[i_edge[i][1]][i_link[i][1]][0];
				}
			}
		}
		//---------- LDPC component ---------- 
		for (int inter_it = 1; inter_it <= inter_it_max; inter_it++)
		{
			//---------- check-to-variable messages ----------
			for (int i = 0; i < h_row; i++)
			{
				for (int j = 0; j < row_w[i]; j++)
				{
					bool flag = true;
					for (int k = 0; k < row_w[i]; k++)
					{
						if (k != j)
						{
							if (flag)
							{
								memcpy(p_c2v[i][row_edge[i][j]], p_v2c[row_edge[i][k]][i], sizeof(double)*NUM_LEVEL);
								flag = false;
							}
							else
							{
								double temp[NUM_LEVEL] = { 0 };
								memcpy(temp, p_c2v[i][row_edge[i][j]], sizeof(double)*NUM_LEVEL);
								if (NUM_LEVEL == 4)
								{
									p_c2v[i][row_edge[i][j]][0] = temp[0] * p_v2c[row_edge[i][k]][i][0] + temp[1] * p_v2c[row_edge[i][k]][i][1] + temp[2] * p_v2c[row_edge[i][k]][i][2] + temp[3] * p_v2c[row_edge[i][k]][i][3];
									p_c2v[i][row_edge[i][j]][1] = temp[0] * p_v2c[row_edge[i][k]][i][1] + temp[1] * p_v2c[row_edge[i][k]][i][0] + temp[2] * p_v2c[row_edge[i][k]][i][3] + temp[3] * p_v2c[row_edge[i][k]][i][2];
									p_c2v[i][row_edge[i][j]][2] = temp[0] * p_v2c[row_edge[i][k]][i][2] + temp[1] * p_v2c[row_edge[i][k]][i][3] + temp[2] * p_v2c[row_edge[i][k]][i][0] + temp[3] * p_v2c[row_edge[i][k]][i][1];
									p_c2v[i][row_edge[i][j]][3] = temp[0] * p_v2c[row_edge[i][k]][i][3] + temp[1] * p_v2c[row_edge[i][k]][i][2] + temp[2] * p_v2c[row_edge[i][k]][i][1] + temp[3] * p_v2c[row_edge[i][k]][i][0];
								}
								else if (NUM_LEVEL == 8)
								{
									p_c2v[i][row_edge[i][j]][0] = temp[0] * p_v2c[row_edge[i][k]][i][0] + temp[1] * p_v2c[row_edge[i][k]][i][1] + temp[2] * p_v2c[row_edge[i][k]][i][2] + temp[3] * p_v2c[row_edge[i][k]][i][3] + temp[4] * p_v2c[row_edge[i][k]][i][4] + temp[5] * p_v2c[row_edge[i][k]][i][5] + temp[6] * p_v2c[row_edge[i][k]][i][6] + temp[7] * p_v2c[row_edge[i][k]][i][7];
									p_c2v[i][row_edge[i][j]][1] = temp[0] * p_v2c[row_edge[i][k]][i][1] + temp[1] * p_v2c[row_edge[i][k]][i][0] + temp[2] * p_v2c[row_edge[i][k]][i][3] + temp[3] * p_v2c[row_edge[i][k]][i][2] + temp[4] * p_v2c[row_edge[i][k]][i][5] + temp[5] * p_v2c[row_edge[i][k]][i][4] + temp[6] * p_v2c[row_edge[i][k]][i][7] + temp[7] * p_v2c[row_edge[i][k]][i][6];
									p_c2v[i][row_edge[i][j]][2] = temp[0] * p_v2c[row_edge[i][k]][i][2] + temp[1] * p_v2c[row_edge[i][k]][i][3] + temp[2] * p_v2c[row_edge[i][k]][i][0] + temp[3] * p_v2c[row_edge[i][k]][i][1] + temp[4] * p_v2c[row_edge[i][k]][i][6] + temp[5] * p_v2c[row_edge[i][k]][i][7] + temp[6] * p_v2c[row_edge[i][k]][i][4] + temp[7] * p_v2c[row_edge[i][k]][i][5];
									p_c2v[i][row_edge[i][j]][3] = temp[0] * p_v2c[row_edge[i][k]][i][3] + temp[1] * p_v2c[row_edge[i][k]][i][2] + temp[2] * p_v2c[row_edge[i][k]][i][1] + temp[3] * p_v2c[row_edge[i][k]][i][0] + temp[4] * p_v2c[row_edge[i][k]][i][7] + temp[5] * p_v2c[row_edge[i][k]][i][6] + temp[6] * p_v2c[row_edge[i][k]][i][5] + temp[7] * p_v2c[row_edge[i][k]][i][4];
									p_c2v[i][row_edge[i][j]][4] = temp[0] * p_v2c[row_edge[i][k]][i][4] + temp[1] * p_v2c[row_edge[i][k]][i][5] + temp[2] * p_v2c[row_edge[i][k]][i][6] + temp[3] * p_v2c[row_edge[i][k]][i][7] + temp[4] * p_v2c[row_edge[i][k]][i][0] + temp[5] * p_v2c[row_edge[i][k]][i][1] + temp[6] * p_v2c[row_edge[i][k]][i][2] + temp[7] * p_v2c[row_edge[i][k]][i][3];
									p_c2v[i][row_edge[i][j]][5] = temp[0] * p_v2c[row_edge[i][k]][i][5] + temp[1] * p_v2c[row_edge[i][k]][i][4] + temp[2] * p_v2c[row_edge[i][k]][i][7] + temp[3] * p_v2c[row_edge[i][k]][i][6] + temp[4] * p_v2c[row_edge[i][k]][i][1] + temp[5] * p_v2c[row_edge[i][k]][i][0] + temp[6] * p_v2c[row_edge[i][k]][i][3] + temp[7] * p_v2c[row_edge[i][k]][i][2];
									p_c2v[i][row_edge[i][j]][6] = temp[0] * p_v2c[row_edge[i][k]][i][6] + temp[1] * p_v2c[row_edge[i][k]][i][7] + temp[2] * p_v2c[row_edge[i][k]][i][4] + temp[3] * p_v2c[row_edge[i][k]][i][5] + temp[4] * p_v2c[row_edge[i][k]][i][2] + temp[5] * p_v2c[row_edge[i][k]][i][3] + temp[6] * p_v2c[row_edge[i][k]][i][0] + temp[7] * p_v2c[row_edge[i][k]][i][1];
									p_c2v[i][row_edge[i][j]][7] = temp[0] * p_v2c[row_edge[i][k]][i][7] + temp[1] * p_v2c[row_edge[i][k]][i][6] + temp[2] * p_v2c[row_edge[i][k]][i][5] + temp[3] * p_v2c[row_edge[i][k]][i][4] + temp[4] * p_v2c[row_edge[i][k]][i][3] + temp[5] * p_v2c[row_edge[i][k]][i][2] + temp[6] * p_v2c[row_edge[i][k]][i][1] + temp[7] * p_v2c[row_edge[i][k]][i][0];
								}
								for (int z = 0; z < NUM_LEVEL; z++)
								{
									if (p_c2v[i][row_edge[i][j]][z] < NUMERIC_LIMIT) p_c2v[i][j][z] = NUMERIC_LIMIT;
								}
							}
						}
					}
				}
			}
			//---------- variable-to-check messages ----------
			for (int i = 0; i < h_col; i++)
			{
				for (int j = 0; j < col_w[i]; j++)
				{
					memcpy(p_v2c[i][col_edge[i][j]], output_p[i], sizeof(double)*NUM_LEVEL);
					for (int k = 0; k < col_w[i]; k++)
					{
						if (k != j)
						{
							for (int z = 0; z < NUM_LEVEL; z++)
							{
								p_v2c[i][col_edge[i][j]][z] *= p_c2v[col_edge[i][k]][i][z];
							}
							double temp = 0;
							for (int z = 0; z < NUM_LEVEL; z++)
							{
								if (p_v2c[i][col_edge[i][j]][z] < NUMERIC_LIMIT) p_v2c[i][col_edge[i][j]][z] = NUMERIC_LIMIT;
								temp += p_v2c[i][col_edge[i][j]][z];
							}
							for (int z = 0; z < NUM_LEVEL; z++)
							{
								p_v2c[i][col_edge[i][j]][z] /= temp;
								if (p_v2c[i][col_edge[i][j]][z] < NUMERIC_LIMIT) p_v2c[i][col_edge[i][j]][z] = NUMERIC_LIMIT;
							}
						}
					}
				}
				for (int j = 0; j < NUM_LEVEL; j++)
				{
					input_p[i][j] = 1;
				}
				for (int j = 0; j < col_w[i]; j++)
				{
					for (int z = 0; z < NUM_LEVEL; z++)
					{
						input_p[i][z] *= p_c2v[col_edge[i][j]][i][z];
					}
					double temp = 0;
					for (int z = 0; z < NUM_LEVEL; z++)
					{
						if (input_p[i][z] < NUMERIC_LIMIT) input_p[i][z] = NUMERIC_LIMIT;
						temp += input_p[i][z];
					}
					for (int z = 0; z < NUM_LEVEL; z++)
					{
						input_p[i][z] /= temp;
						if (input_p[i][z] < NUMERIC_LIMIT) input_p[i][z] = NUMERIC_LIMIT;
					}
				}
			}
		}
	}
	//---------- finalization ----------
	for (int i = 0; i < h_col; i++)
	{
		memcpy(app[i], output_p[i], sizeof(double)*NUM_LEVEL);
		for (int j = 0; j < col_w[i]; j++)
		{
			for (int z = 0; z < NUM_LEVEL; z++)
			{
				app[i][z] *= p_c2v[col_edge[i][j]][i][z];
			}
			double temp = 0;
			for (int z = 0; z < NUM_LEVEL; z++)
			{
				if (app[i][z] < NUMERIC_LIMIT) app[i][z] = NUMERIC_LIMIT;
				temp += app[i][z];
			}
			for (int z = 0; z < NUM_LEVEL; z++)
			{
				app[i][z] /= temp;
				if (app[i][z] < NUMERIC_LIMIT) app[i][z] = NUMERIC_LIMIT;
			}
		}
	}
}