#include<iostream>
#include<vector>
#include"lt.h"
#include"degree_profile.h"
#include<algorithm>
#include<fstream>
#define e 2.718281828

#define ROUND(x) int(x+0.5)
using namespace std;

void create_G(vector<vector<int>>&, int, int, vector<int>&, vector<int>&, vector<int>&, vector<int>&,int &);


LT::LT(int inter_code_len, int max_code_len)
{
	this->inter_code_len = inter_code_len;
	this->max_code_len = max_code_len;
	vector<int>inter_code(inter_code_len, 0);
	vector<int>max_code(max_code_len, 0);
	G.resize(inter_code_len, max_code);

	column_bits.resize(max_code_len, 0);					//     bits/column
	row_bits.resize(inter_code_len, 0);						//     bits/row
	column_address.resize(max_code_len * 20);				//	   the address in column.According to different H, the size will be different.
	row_address.resize(max_code_len * 20);					//     the address in row

	create_G(G, inter_code_len, max_code_len, column_bits, row_bits, column_address, row_address,total_number);

}

vector<int> LT::encoder(vector<int> tx_ldpc, int recursive)
{
	if (SYSTEMATIC)
	{
		if (recursive == 0)
		{
			return tx_ldpc;
		}
		else
		{
			vector<int> added_codeword(IR_LEN);

			int reg = 0, reg2 = 0;
		/*	for (int i = 0; i < IR_LEN; i++)
			{
				reg = 0;
				//cout << inter_code_len + (recursive - 1)*IR_LEN + i<<",";
				for (int j = 0; j < inter_code_len; j++)
				{
					reg = reg + G[j][inter_code_len + (recursive - 1)*IR_LEN + i] * tx_ldpc[j];
				}
				added_codeword[i] = reg % 2;
			}*/
			for (int i = 0; i < inter_code_len + (recursive - 1)*IR_LEN; i++)
			{
				reg2 += column_bits[i];
			}
			for (int i = inter_code_len + (recursive - 1)*IR_LEN; i < inter_code_len + recursive*IR_LEN; i++)
			{
				reg = 0;
				for (int j=0; j < column_bits[i] ; j++)
				{
					reg += tx_ldpc[column_address[reg2 + j]-1];
				}
				added_codeword[i- inter_code_len - (recursive - 1)*IR_LEN] = reg % 2;
				reg2 += column_bits[i];
			}

			return added_codeword;

			added_codeword.clear();
		}
	}
	else
	{
		if (recursive == 0)
		{
			vector<int> origin_codeword(INTER_CODE_LEN);

			int reg = 0, reg2 = 0;
			
			for (int i = 0; i < inter_code_len; i++)
			{
				reg = 0;
				for (int j = 0; j < column_bits[i]; j++)
				{
					reg += tx_ldpc[column_address[reg2 + j] - 1];
				}
				origin_codeword[i] = reg % 2;
				reg2 += column_bits[i];
			}

			return origin_codeword;

			origin_codeword.clear();
		}
		else
		{
			vector<int> added_codeword(IR_LEN);

			int reg = 0, reg2 = 0;
		
			for (int i = 0; i < inter_code_len + (recursive - 1)*IR_LEN; i++)
			{
				reg2 += column_bits[i];
			}
			for (int i = inter_code_len + (recursive - 1)*IR_LEN; i < inter_code_len + recursive * IR_LEN; i++)
			{
				reg = 0;
				for (int j = 0; j < column_bits[i]; j++)
				{
					reg += tx_ldpc[column_address[reg2 + j] - 1];
				}
				added_codeword[i - inter_code_len - (recursive - 1)*IR_LEN] = reg % 2;
				reg2 += column_bits[i];
			}

			return added_codeword;

			added_codeword.clear();
		}
	}
}

vector<double> LT::decoder(double snr, vector<double>rx, int recursive)
{
	
	if (recursive < 1)
	{
		vector<double> r(inter_code_len);
		for (int i = 0; i < inter_code_len; i++)
			r[i] = 4 * rx[i] *pow(10, snr / 10);

		return r;
		r.clear();
	}
	else
	{
		vector<double> r(inter_code_len + recursive * IR_LEN);
		vector<int> z(inter_code_len + recursive * IR_LEN);
		vector<vector<double>> M(inter_code_len, r);
		vector<vector<double>> E(inter_code_len, r);
		vector<double> L(inter_code_len, 0);
		vector<double> llr_out(inter_code_len);
		bool triger = true;
		int iteration = 0;
		double reg = 0;
		int reg2 = 0;

		

		for (int i = 0; i < inter_code_len + recursive * IR_LEN; i++)
		{
			//r[i] = 4 * rx[i] * 0.5*pow(10, snr / 10);
			r[i] = 4 * rx[i] *pow(10, snr / 10);
		}
		
		while (iteration < 20)
		{
			if (iteration == 0)
			{
				//for (int j = 0; j < inter_code_len + recursive * IR_LEN; j++)
				for (int j = 0; j < inter_code_len; j++)
				{
					for (int i = 0; i < column_bits[j]; i++)
					{
						E[column_address[i + reg2] - 1][j] =r[j];
					//	E[column_address[i + reg2] - 1][j] = log((1 + tanh(r[j] / 2)) / (1 - tanh(r[j] / 2)));
						if (E[column_address[i + reg2] - 1][j] > 30)
						{
							E[i][row_address[j + reg2] - 1] = 40;
						}
						else if (E[column_address[i + reg2] - 1][j] < -30)
						{
							E[i][row_address[j + reg2] - 1] = -40;
						}
					}
					reg2 += column_bits[j];
				}
			}
			else
			{
				reg = 1;
				for (int i = 0; i < inter_code_len + recursive * IR_LEN; i++)
				{
					for (int j = 0; j < column_bits[i]; j++)
					{
						for (int k = 0; k < column_bits[i]; k++)
						{
							if (k != j)
								reg *= tanh(M[column_address[k + reg2] - 1][i] / 2);
						}
						reg *= tanh(r[i] / 2);

						E[column_address[j + reg2] - 1][i] = log((1 + reg) / (1 - reg));
						if (E[column_address[j + reg2] - 1][i] > 30)
						{
							E[column_address[j + reg2] - 1][i] = 40;
						}
						else if (E[column_address[j + reg2] - 1][i] < -30)
						{
							E[column_address[j + reg2] - 1][i] = -40;
						}
						reg = 1;
					}
					reg2 += column_bits[i];
				}

			}

			reg = 0;
			reg2 = 0;
			for (int i = 0; i < inter_code_len; i++)
			{
				for (int j = 0; j < row_bits[i] && row_address[reg2 + j] < inter_code_len + recursive * IR_LEN; j++)
				{
					for (int k = 0; k < row_bits[i] && row_address[reg2 + k] < inter_code_len + recursive * IR_LEN; k++)
					{
						if (k != j)
						{
							reg += E[i][row_address[k + reg2] - 1];
						}
					}
					M[i][row_address[j + reg2] - 1] = reg;
					reg = 0;
				}
				reg2 += row_bits[i];
			}
			reg = 0;
			reg2 = 0;

			iteration++;
		}

		for (int i = 0; i < inter_code_len; i++)
		{
			for (int j = 0; j < row_bits[i] && row_address[reg2 + j]< inter_code_len + recursive * IR_LEN; j++)
			{
				L[i] += E[i][row_address[reg2 + j] - 1];
			}
			reg2 += row_bits[i];
		}

		reg = 0;
		reg2 = 0;

		

		return L;
		z.clear();
		r.clear();
		M.clear();
		E.clear();
		L.clear();

	}
	
}

void create_G(vector<vector<int>>&G, int row_G, int column_G, vector<int>&column_bits, vector<int>&row_bits, vector<int>&column_address, vector<int>&row_address,int & total_number)
{
	vector<int>o_edge(column_G);

	if (SYSTEMATIC)
	{
		//--select the degree of output check nodes--
		for (int i = 0; i < row_G; i++)
			o_edge[i] = 1;

		int reg = 0;
		int new_step[STEP_SIZE + 1];

		if (STEP_SIZE == 3)
		{
			new_step[0] = row_G;
			new_step[1] = step[0];
			new_step[2] = step[1];
			new_step[3] = column_G;
		}
		else
		{
			new_step[0] = row_G;
			new_step[1] = column_G;
		}

		for (int l = 0; l < STEP_SIZE; l++)
		{
			for (int i = 0; i < DEG_NUM; i++)
			{
				int j = 0;
				while (j < ROUND(deg_pdf[l][i] * (new_step[l + 1] - new_step[l])))
				{
					o_edge[reg + row_G] = deg[i];
					reg++;
					j++;
				}
			}
		}

		//--interleaver of output check nodes
		reg = 0;
		int reg2 = 0;
		for (int l = 0; l < STEP_SIZE; l++)
		{
			for (int i = new_step[l]; i < new_step[l + 1]; i++)
			{
				reg2 = new_step[l] + rand() % (new_step[l + 1] - new_step[l]);
				reg = o_edge[i];
				o_edge[i] = o_edge[reg2];
				o_edge[reg2] = reg;
			}
		}	
		//--create G--
		reg = 0, reg2 = 0;
		for (int i = 0; i < row_G; i++)
			G[i][i] = 1;

		vector<int> non_repeat(row_G);
		int t = 0;
		generate(non_repeat.begin(), non_repeat.end(), [&] {return t++; });
		for (int i = row_G; i < column_G; i++)
		{
			reg = 0, reg2 = 0;
			for (int k = 0; k < o_edge[i]; k++)
			{
				reg2 = rand() % row_G;
				reg = non_repeat[k];
				non_repeat[k] = non_repeat[reg2];
				non_repeat[reg2] = reg;
			}
			for (int j = 0; j < o_edge[i]; j++)
			{
				G[non_repeat[j]][i] = 1;
			}
		}
	}
	else
	{
		//--select the degree of output check nodes--
		int reg = 0;
		int new_step[4];
		
		new_step[0] = 0;
		new_step[1] = step[0];
		new_step[2] = step[1];
		new_step[3] = column_G;

		for (int l = 0; l < STEP_SIZE; l++)
		{
			for (int i = 0; i < DEG_NUM; i++)
			{
				int j = 0;
				while (j < ROUND(deg_pdf[l][i] * (new_step[l + 1] - new_step[l])))
				{
					o_edge[reg] = deg[i];
					reg++;
					j++;
				}
			}
		}

		//--interleaver of output check nodes
		reg = 0;
		int reg2 = 0;
		for (int l = 0; l < STEP_SIZE; l++)
		{
			for (int i = new_step[l]; i < new_step[l + 1]; i++)
			{
				reg2 = new_step[l] + rand() % (new_step[l + 1] - new_step[l]);
				reg = o_edge[i];
				o_edge[i] = o_edge[reg2];
				o_edge[reg2] = reg;
			}
		}
		//--create G--
		reg = 0, reg2 = 0;

		vector<int> non_repeat(column_G);
		int t = 0;
		generate(non_repeat.begin(), non_repeat.end(), [&] {return t++; });
		for (int i = 0; i < column_G; i++)
		{
			reg = 0, reg2 = 0;
			for (int k = 0; k < o_edge[i]; k++)
			{
				reg2 = rand() % row_G;
				reg = non_repeat[k];
				non_repeat[k] = non_repeat[reg2];
				non_repeat[reg2] = reg;
			}
			for (int j = 0; j < o_edge[i]; j++)
			{
				G[non_repeat[j]][i] = 1;
			}
		}
	}

	o_edge.clear();
	//--input:G ,output:address bits
	int reg = 0, reg2 = 0;
	for (int i = 0; i < column_G; i++)
	{
		for (int j = 0; j < row_G; j++)
		{
			if (G[j][i] == 1)
			{
				column_address[reg2] = j + 1;
				reg++;
				reg2++;
			}
		}
		column_bits[i] = reg;
		reg = 0;
	}

	reg2 = 0;
	for (int i = 0; i < row_G; i++)
	{
		for (int j = 0; j < column_G; j++)
		{
			if (G[i][j] == 1)
			{
				row_address[reg2] = j + 1;
				reg++;
				reg2++;
			}
		}
		row_bits[i] = reg;
		reg = 0;
	}

	total_number = reg2;
	//--ofile G
	ofstream ofile("LT_G.txt", ios::out);

	if (ofile.is_open())
	{
		for (auto val : column_bits)
			ofile << val << " ";
		ofile << endl << endl;
		for (auto val : row_bits)
			ofile << val << " ";
		ofile << endl << endl;

		for (int i = 0; i < reg2; i++)
			ofile << column_address[i] << " ";
		ofile << endl << endl;
		for (int i = 0; i < reg2; i++)
			ofile << row_address[i] << " ";
		ofile << endl << endl;
	}
	else
		cout << "¼g¤J¥¢±Ñ";
	ofile.close();
}

