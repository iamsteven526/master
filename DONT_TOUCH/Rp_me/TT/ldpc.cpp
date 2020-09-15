#include <iostream>
#include "ldpc.h"
#include <fstream>
#include <algorithm>


using namespace std;


void print_matrix(vector<vector<int>>&,int,int,char);
void create_G(vector<vector<int>>&, vector<vector<int>>&, int, int,vector<int>&,vector<int>&, vector<int>&, vector<int>&,int &);
void add(int, int, int, vector<vector<int>>&);
void change(int, int, int, vector<vector<int>>&);

LDPC::LDPC(int column, int row)
{
	message.resize(column-row, 0);
	codeword.resize(column, 0);
	G.resize(column-row, codeword);
	H.resize(row, codeword);
	this->column = column;
	this->row = row;

	column_address.resize(column*3);	//	   the address in column.According to different H, the size will be different.
	row_address.resize(column*3);		//     the address in row
	column_bits.resize(column, 0);		//     bits/column
	row_bits.resize(row, 0);			//     bits/row

	create_G(G, H, row, column,row_address,column_address,row_bits,column_bits,total_number);
	//print_matrix(H,row, column, 'H');
	//print_matrix(G,column-row,column,'G');
}

vector<int> LDPC::encoder(vector<int> random_message)
{
//	vector<int> random_message(column-row);
//	generate(random_message.begin(), random_message.end(), [] {return rand()%2; });
//	generate(random_message.begin(), random_message.end(), [] {return 1; });

	int reg;
	for (int i = 0; i < column; i++)
	{
		reg = 0;
		for (int j = 0; j < column-row; j++)
		{
			reg = reg + G[j][i] * random_message[j];
		}
		codeword[i] = reg % 2;
	}
	random_message.clear();

	return codeword;
}

vector<double> LDPC::decoder(double snr,vector<double>&rx,bool&out_triger)
{
	vector<double> r(column,0);
	vector<int> z(column,0);
	vector<vector<double>> M(row, r);
	vector<vector<double>> E(row, r);
	vector<double> L(column);
	bool triger = true;
	int iteration = 0;
	double reg = 0;
	int reg2 = 0;
	int reg3 = 0;
	
	for (int i = 0; i < column; i++)
	{
		r[i] = rx[i];
	}

	while (triger == true && iteration < 10)
	{
		triger = false;
		if (iteration == 0)
		{
			for (int j = 0; j < column; j++)
			{
				for (int i = 0; i < column_bits[j]; i++)
				{
					M[column_address[i+reg2] - 1][j] = r[j];
				}
				reg2 += column_bits[j];
			}
		}
		else
		{
			for (int i = 0; i < column; i++)
			{
				for (int j = 0; j < column_bits[i]; j++)
				{
					for (int k = 0; k < column_bits[i]; k++)
					{
						if(k!=j)
							reg += E[column_address[k + reg2] - 1][i];
					}
					M[column_address[j + reg2] - 1][i] = r[i] + reg;
					reg = 0;
				}
				reg2 += column_bits[i];
			}
		}

		reg = 1;
		reg2 = 0;
		for (int i = 0; i < row; i++)
		{
			for (int j = 0; j < row_bits[i]; j++)
			{
				for (int k = 0; k < row_bits[i]; k++)
				{
					if(k!=j)
						reg *= tanh(M[i][row_address[k + reg2]-1] / 2);
				}
				E[i][row_address[j + reg2]-1] = (log((1 + reg) / (1 - reg)));
				reg = 1;
			}
			reg2 += row_bits[i];
		}
	
		reg = 0;
		reg2 = 0;

		for (int i = 0; i < column; i++)
		{
			for (int j = 0; j < column_bits[i]; j++)
			{
				reg += E[column_address[j + reg2]-1][i];
			}
			L[i] = r[i] + reg;
			reg = 0;
			reg2 += column_bits[i];

			if (L[i] < 0)
				z[i] = 1;
			else
				z[i] = 0;
		}
		
		reg2 = 0;
	
		for (int i = 0; i < row; i++)
		{
			for (int j = 0; j < row_bits[i]; j++)
			{
				reg3=reg3+z[row_address[reg2]-1];
				reg2++;
			}
			if (reg3 % 2 != 0)
			{
				triger = true;
			}
			reg3 = 0;
		}
		reg2 = 0;
		reg3 = 0;
		if (triger == false)
			out_triger = false;

		iteration++ ;
	}
	
	return L;
	z.clear();
	r.clear();
	M.clear();
	E.clear();
	L.clear();
	
}

void create_G(vector<vector<int>>&G, vector<vector<int>>&H, int row, int column,vector<int>&row_address,vector<int>&column_address,vector<int>&row_bits,vector<int>&column_bits,int & total_number)
{
	///////////////////////////read
	int in;
	
	ifstream fin("H_2000_100.txt");
	if (!fin)
	{
		cout << "can't read";
		system("pause");
	}

	for (int i = 0; i < column; i++)
	{
		fin >> in;
		column_bits[i] = in;
	}

	for (int i = 0; i < row; i++)
	{
		fin >> in;
		row_bits[i] = in;
	}
	int reg = 0;
	for (int i = 0; i < column; i++)
	{
		for (int j = 0; j < column_bits[i]; j++)
		{
			fin >> in;
			H[in - 1][i] = 1;
			column_address[reg + j] = in;
		}
		reg += column_bits[i];
	}

	reg = 0;
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < row_bits[i]; j++)
		{
			fin >> in;
			row_address[reg + j]=in;
		}
		reg += row_bits[i];
	}

	total_number = reg;
	fin.close();



	//////////////////////////gaussian
	for (int i = 0; i < min(column, row); i++)
	{
		if (H[i][i] == 0)
		{
			int k = 0;
			for (k = i + 1; k < row; k++)
			{
				if (H[k][i] == 1)
				{
					change(k, i, column, H);
					break;
				}
			}
			for (k = i + 1; k < row; k++)
			{
				if (H[k][i] == 1)
				{
					add(i, k, column, H);
				}
			}
		}
		else
		{
			for (int k = i + 1; k < row; k++)
			{
				if (H[k][i] == 1)
				{
					add(i, k, column, H);
				}
			}
		}
	}

	for (int i = min(column, row) - 1; i >0; i--)
	{
		if (H[i][i] == 1)
		{
			for (int k = i - 1; k >= 0; k--)
			{
				if (H[k][i] == 1)
				{
					add(i, k,column,H);
				}
			}
		}
	}

	///////////////////////////transpose
	for (int i = 0; i <column-row; i++)
	{
		for (int j = 0; j < row; j++)
			G[i][j] = H[j][row+i];
	}

	for (int i = row; i < column; i++)
		G[i - row][i] = 1;
	
}

void change(int a, int b,int column, vector<vector<int>>&H)
{
	int reg = 0;

	for (int i = 0; i < column; i++)
	{
		reg = H[a][i];
		H[a][i] = H[b][i];
		H[b][i] = reg;
	}
}

void add(int a, int b, int column, vector<vector<int>>&H)
{
	for (int i = 0; i < column; i++)
	{
		H[b][i] = (H[a][i] + H[b][i]) % 2;
	}
}

void print_matrix(vector<vector<int>>&matrix,int row,int column,char matrix_name)
{
	ofstream ofile("HH.txt", ios::out);

	if (ofile.is_open())
	{
		ofile << matrix_name << endl;
		for (int i = 0; i < row; i++)
		{
			for (int j = 0; j < column; j++)
			{
				ofile << matrix[i][j];
			}
			ofile << endl;
		}
	}
	else
		cout << "¼g¤J¥¢±Ñ";

	ofile.close();
}