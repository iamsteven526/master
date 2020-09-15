#include <vector>
#include <iostream>
#include "raptor.h"
#include <fstream>
#include "parameter.h"

using namespace std;


Raptor::Raptor(LDPC&ldpc, LT&lt,int data_len,int inter_code_len,int max_code_len)
{
	this->data_len = data_len;
	this->inter_code_len = inter_code_len;
	this->max_code_len = max_code_len;

	// the mean of ldpc's column_bits and column_address is different to the lt's

	column_bits.reserve(ldpc.row_bits.size() + lt.column_bits.size());
	column_bits.insert(column_bits.end(), ldpc.row_bits.begin(), ldpc.row_bits.end());
	column_bits.insert(column_bits.end(), lt.column_bits.begin(), lt.column_bits.end());
	column_address.reserve(ldpc.total_number + lt.total_number);
	column_address.insert(column_address.end(), ldpc.row_address.begin(), ldpc.row_address.begin() + ldpc.total_number);
	column_address.insert(column_address.end(), lt.column_address.begin(), lt.column_address.begin() + lt.total_number);

	row_bits.resize(inter_code_len);
	row_address.resize(ldpc.total_number + lt.total_number);


	for (int i = 0; i < inter_code_len; i++)
	{
		row_bits[i] = ldpc.column_bits[i] + lt.row_bits[i];
	}
	

	int reg=0,reg2 = 0,reg3=0;
	for (int i = 0; i < inter_code_len; i++)
	{
		for (int j = 0; j < ldpc.column_bits[i]; j++)
		{
			row_address[reg3 + j] = ldpc.column_address[j + reg];
		}
		reg += ldpc.column_bits[i];
		reg3 += ldpc.column_bits[i];

		for (int j = 0; j < lt.row_bits[i]; j++)
		{
			row_address[reg3 + j] = lt.row_address[j + reg2]+inter_code_len-data_len;
		}
		reg2 += lt.row_bits[i];
		reg3+=lt.row_bits[i];

	}
	reg = 0;
	reg2 = 0;
	reg3 = 0;


	//--ofile joint matrix
	ofstream ofile("Joint_matrix.txt", ios::out);

	if (ofile.is_open())
	{
		for (auto val : column_bits)
			ofile << val << " ";
		ofile << endl << endl;
		for (auto val : row_bits)
			ofile << val << " ";
		ofile << endl << endl;

		for (int i = 0; i < ldpc.total_number + lt.total_number; i++)
			ofile << column_address[i] << " ";
		ofile << endl << endl;
		for (int i = 0; i < ldpc.total_number + lt.total_number; i++)
			ofile << row_address[i] << " ";
		ofile << endl << endl;
	}
	else
		cout << "¼g¤J¥¢±Ñ";
	ofile.close();

}

vector<double> Raptor::joint_decoder(int recursive, double snr, vector<double> rx, bool&out_triger,LDPC&ldpc)
{
	if (recursive < 1)
	{
		vector<double> r(INTER_CODE_LEN);
		for (int i = 0; i < INTER_CODE_LEN; i++)
			r[i] = 4 * rx[i] * 0.5*pow(10, snr / 10);

		return  ldpc.decoder(snr, r, out_triger);
		r.clear();
	}
	else
	{
		vector<double> r(INTER_CODE_LEN + recursive * IR_LEN + (INTER_CODE_LEN - DATA_LEN), 0); // we define the previous codes are ldpc check nodes 
		vector<int> z(INTER_CODE_LEN);
		vector<vector<double>> M(INTER_CODE_LEN, r);
		vector<vector<double>> E(INTER_CODE_LEN, r);
		vector<double> L(INTER_CODE_LEN, 0);
		
		bool triger = true;
		int iteration = 0;
		double reg = 0;
		int reg2 = 0;

		for (int i = 0; i < INTER_CODE_LEN + recursive * IR_LEN; i++)
		{
			r[i + (INTER_CODE_LEN - DATA_LEN)] = 4 * rx[i] *pow(10, snr / 10);
		}
		
	//	for (auto val : r)
	//		cout << val<<" ";

	//	system("pause");
	//	cout << log((1 + tanh(r[2026] / 2)) / (1 - tanh(r[2026] / 2))) / log(e) << endl;

		while (iteration < 30 && triger==true)
		{
			triger = false;

			if (iteration == 0)
			{
				//reg2 must begin from (INTER_CODE_LEN - DATA_LEN)-th column
				//not from zero-th column
				for (int j = 0; j < INTER_CODE_LEN - DATA_LEN; j++)
					reg2 += column_bits[j];

				for (int j = INTER_CODE_LEN - DATA_LEN; j < INTER_CODE_LEN + recursive * IR_LEN + (INTER_CODE_LEN - DATA_LEN); j++)
				{
					for (int i = 0; i < column_bits[j]; i++)
					{
						E[column_address[i + reg2] - 1][j] = log((1 + tanh(r[j] / 2)) / (1 - tanh(r[j] / 2)));

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
				reg2 = 0;



		//		cout << E[0][2026];
		//		system("pause");
			}
			else
			{
				reg = 1;
				for (int i = 0; i < INTER_CODE_LEN + recursive * IR_LEN + (INTER_CODE_LEN - DATA_LEN); i++)
				{
					for (int j = 0; j < column_bits[i]; j++)
					{
						for (int k = 0; k < column_bits[i]; k++)
						{
							if (k != j)
							{
								/*if (tanh(M[column_address[k + reg2] - 1][i] / 2) == 0)
								{
									cout << "error!";
									cout << k + reg2 <<" " <<i<<endl;
								}*/
								reg *= tanh(M[column_address[k + reg2] - 1][i] / 2);
							}
						}
						if(i>=INTER_CODE_LEN-DATA_LEN)
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
			for (int i = 0; i < INTER_CODE_LEN; i++)
			{
				for (int j = 0; j < row_bits[i] && row_address[reg2 + j] < INTER_CODE_LEN + recursive * IR_LEN + (INTER_CODE_LEN - DATA_LEN); j++)
				{
					for (int k = 0; k < row_bits[i] && row_address[reg2 + k] < INTER_CODE_LEN + recursive * IR_LEN + (INTER_CODE_LEN - DATA_LEN); k++)
					{
						if (k != j)
						{
							reg += E[i][row_address[k + reg2] - 1];
						//	cout << row_address[k + reg2] - 1 << " ";
						}
					}
					M[i][row_address[j + reg2] - 1] = reg;
					//if (reg == 0)
					//{
					//	cout << i<<" "<< row_address[j + reg2] - 1<<endl;
					//}
					reg = 0;
				}
				reg2 += row_bits[i];
			}
			reg = 0;
			reg2 = 0;

			for (int i = 0; i < INTER_CODE_LEN; i++)
			{
				for (int j = 0; j < row_bits[i] && row_address[reg2 + j]< INTER_CODE_LEN + recursive * IR_LEN + (INTER_CODE_LEN - DATA_LEN); j++)
				{
					L[i] += E[i][row_address[reg2 + j] - 1];
				//	cout << E[i][row_address[reg2 + j] - 1] << "+";
				}
				//cout << endl;
				reg2 += row_bits[i];

				if (L[i] < 0)
					z[i] = 1;
				else
					z[i] = 0;
			}

			reg = 0;
			reg2 = 0;
			int reg3 = 0;

			// early termiation
			for (int i = 0; i < INTER_CODE_LEN-DATA_LEN; i++) //"INTER_CODE_LEN-DATA_LEN" is the row of LDPC H 
			{
				for (int j = 0; j < ldpc.row_bits[i]; j++)
				{
					reg3 = reg3 + z[ldpc.row_address[reg2] - 1];
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
			reg = 0;
			if (triger == false)
				out_triger = false;

			//cout << triger << " ";
			//cout <<"iteration"<< iteration << endl;
			iteration++;
		}


		return L;

		z.clear();
		r.clear();
		M.clear();
		E.clear();
		L.clear();
	}
}

vector<double> Raptor::tandem_decoder(double snr, vector<double> rx, int recursive, bool & out_triger, LDPC& ldpc, LT& lt)
{
	vector<double>rx_lt = lt.decoder(snr, rx, recursive);
	return ldpc.decoder(snr, rx_lt, out_triger);
}
