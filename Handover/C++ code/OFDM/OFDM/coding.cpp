#include "channel_coding.h"

#define min(x,y) ((x) < (y) ? (x) : (y));

int **LDPC::CreateGeneratorMatrix(vector<uint16_t> &row_operation_order)
{
	BitMatrix matrix;
	
	int **G_P = matrix.CreateMatrix(infoLength, blockLength - infoLength);
	cout << endl << "[Constructing generator matrix...]" << endl;

	int **H_T = matrix.CreateMatrix(blockLength ,blockLength - infoLength);
	int **H_T_new = matrix.CreateMatrix(blockLength, blockLength - infoLength);;
	int **H_new = matrix.CreateMatrix(blockLength - infoLength, blockLength);


	FILE *Hfile;
	Hfile = fopen(H_FILE, "r");
	for (int row_idx = 0; row_idx < blockLength - infoLength; row_idx++)
	{
		for (int col_idx = 0; col_idx < blockLength; col_idx++)
			fscanf(Hfile, "%d", &H_T[col_idx][row_idx]);
	}fclose(Hfile);

	row_operation_order = matrix.FindIndependentRow(H_T);

	for (int i = 0; i < blockLength; i++)
		memcpy(H_T_new[i], H_T[row_operation_order.at(i)], (blockLength - infoLength) * sizeof(int));   //«eN-K­Óindependent row

	for (int i = 0; i < blockLength - infoLength; i++)
	{
		for (int j = 0; j < blockLength; j++)
			H_new[i][j] = H_T_new[j][i];
	}


	int length = min(blockLength - infoLength, blockLength);

	for (int row_idx = 0; row_idx < length; row_idx++)
	{
		for (int idx = row_idx; idx < (blockLength - infoLength); idx++)
		{
			if (H_new[idx][row_idx] == 1)
			{
				for (int col_idx = 0; col_idx < blockLength; col_idx++)
					swap(H_new[row_idx][col_idx], H_new[idx][col_idx]);
				break;
			}
		}
		for (int idx = 0; idx < (blockLength - infoLength); idx++)
		{
			if (idx != row_idx && H_new[idx][row_idx] == 1) {
				for (int col_idx = 0; col_idx < blockLength; col_idx++)
				{
					H_new[idx][col_idx] = H_new[idx][col_idx] ^ H_new[row_idx][col_idx];
				}
			}
		}
	}


	for (int i = 0; i < infoLength; i++)
	{
		for (int j = 0; j < (blockLength - infoLength); j++)
			G_P[i][j] = H_new[j][i + (blockLength - infoLength)];
	}
	matrix.DeleteMatrix(H_T);
	matrix.DeleteMatrix(H_T_new);
	matrix.DeleteMatrix(H_new);

	cout << "[done.]" << endl << endl;
	return G_P;
}

vector<int> LDPC::SubBlockEncoding(vector<int> subInfoBits)
{
	vector<int> codedBits(blockLength);
	vector<int>	orderCodedBits(blockLength);

	for (int i = 0; i < infoLength; i++)
		codedBits.at(i + (blockLength - infoLength)) = subInfoBits.at(i);

	for (int i = 0; i < (blockLength - infoLength); i++)
	{
		int sum = 0;
		for (int idx = 0; idx < infoLength; idx++)
			sum = sum ^ (subInfoBits.at(idx) & G_P[idx][i]);
		codedBits.at(i) = sum;
	}

	for (int i = 0; i < blockLength; i++)
		orderCodedBits.at(i) = codedBits.at(colOperationOrder.at(i));

	return orderCodedBits;
}

vector<int>  LDPC::SubBlockDecoding(vector<double> receivedLLR, vector<double> &subDecodedLLR, int orderIdx, int &err, bool &check)
{
	vector<double> decodedLdpcLLR(blockLength); subDecodedLLR.resize(decodedLdpcLLR.size());
	vector<int> decodedBits(blockLength), orderedDecodedBits(blockLength);

	double **V2C = new double*[(blockLength - infoLength)];
	double **C2V = new double*[(blockLength - infoLength)];
	for (int i = 0; i < (blockLength - infoLength); i++)
	{
		V2C[i] = new double[blockLength];
		memset(V2C[i], 0, blockLength * sizeof(double));

		C2V[i] = new double[blockLength];
		memset(C2V[i], 0, blockLength * sizeof(double));
	}

	/*  first iteration - initial Channel LLR  */
	for (int col = 0; col < blockLength; col++)
	{
		for (int index = 0; index < checkNodeSet[col].size(); index++)
			V2C[checkNodeSet[col].at(index)][col] = receivedLLR.at(col);
	}

	 int iteration = 0;
	while (check == false && iteration < ITERATION)
	{
		iteration++;

		for (int row = 0; row < (blockLength - infoLength); row++)
		{
			for (int idx = 0; idx < variableNodeSet[row].size(); idx++)
			{
				double tanhvalue = 1.0;
				for (int inner_idx = 0; inner_idx < variableNodeSet[row].size(); inner_idx++)
				{
					if (inner_idx != idx)
						tanhvalue = tanhvalue * tanh(V2C[row][variableNodeSet[row].at(inner_idx)] / 2);
				}
				if (atanh(tanhvalue) > TANH_BOUND)
					C2V[row][variableNodeSet[row].at(idx)] = 2 * TANH_BOUND;
				else if (atanh(tanhvalue) < -TANH_BOUND)
					C2V[row][variableNodeSet[row].at(idx)] = 2 * (-TANH_BOUND);
				else
					C2V[row][variableNodeSet[row].at(idx)] = 2 * atanh(tanhvalue);
			}
		}

		for (int col = 0; col < blockLength; col++)
		{

			for (int index = 0; index < checkNodeSet[col].size(); index++)
			{
				double sum = 0.0;
				for (int inner_idx = 0; inner_idx < checkNodeSet[col].size(); inner_idx++)
				{
					if (inner_idx != index)
						sum = sum + C2V[checkNodeSet[col].at(inner_idx)][col];
				}
				V2C[checkNodeSet[col].at(index)][col] = sum + receivedLLR.at(col);
			}

			double sum = 0.0;
			for (int index = 0; index < checkNodeSet[col].size(); index++)
				sum = sum + C2V[checkNodeSet[col].at(index)][col];

			/*if (abs(receivedLLR.at(col)) > TANH_BOUND)
				decodedLdpcLLR.at(col) = sum + (receivedLLR.at(col) / abs(receivedLLR.at(col)))*(TANH_BOUND + 10);
			else
				decodedLdpcLLR.at(col) = sum + receivedLLR.at(col);*/
			decodedLdpcLLR.at(col) = sum + receivedLLR.at(col);
			

			(decodedLdpcLLR.at(col) >= 0) ? decodedBits.at(col) = 0 : decodedBits.at(col) = 1;
			subDecodedLLR.at(col) = decodedLdpcLLR.at(col);
		}


		vector<int> syndrome(blockLength - infoLength, 0);
		for (int r = 0; r < blockLength; r++)
		{
			if (decodedBits.at(r) == 1)
			{
				for (int c = 0; c < blockLength - infoLength; c++)
					syndrome.at(c) ^= H[c][r];
			}
		}

		int sum = 0;
		for (int i = 0; i < syndrome.size(); i++)
			sum += syndrome.at(i);
		if (sum == 0) check = 1;
	}


	for (int i = 0; i < blockLength; i++)
		orderedDecodedBits.at(colOperationOrder.at(i)) = decodedBits.at(i);

	vector<double> orderLLR(blockLength);
	for (int i = 0; i < blockLength; i++)
		orderLLR.at(colOperationOrder.at(i)) = decodedLdpcLLR.at(i);

	vector<double> inOrderLLR(blockLength);
	for (int i = 0; i < blockLength; i++)
		inOrderLLR.at(colOperationOrder.at(i)) = receivedLLR.at(i);


	for (int i = 0; i < (blockLength - infoLength); i++)
	{
		delete[] C2V[i];
		delete[] V2C[i];
	}
	delete[]C2V; delete[] V2C;
	
	err = 0;


	for (int i = 0; i < infoLength; i++)
	{
		if (orderedDecodedBits.at(i + (blockLength - infoLength)) != infoBits.at(infoLength*orderIdx + i))
		{
			err++;
			//cout <<endl << i <<", output LLR: " << orderLLR.at(i + (blockLength - infoLength))<<", in LLR: "<<receivedLLR.at(i + (blockLength - infoLength));
		}
	}
	/*if (err != 0) 
	{
		system("pause");
		cout << endl;
	}*/

	/*double sum = 0;
	for (auto &j : receivedLLR)
		sum += abs(j);
	cout << endl;
	cout << "ave LLR: " << sum / receivedLLR.size() << ", err: " << err << endl;*/

	
	return decodedBits;
}

vector<uint16_t> LDPC::CreateRow2colOrder(vector<uint16_t> row_operation_order)
{
	int length = row_operation_order.size();
	vector<uint16_t> col_operation_order(length);

	for (int index = 0; index < length; index++)
	{
		for (int i = 0; i < length; i++)
		{
			if (row_operation_order.at(i) == index)
			{
				col_operation_order.at(index) = i;
				break;
			}
		}
	}

	return col_operation_order;
}

vector<int> LDPC::MutlipleSubBlockEncoder()
{
	infoBits.resize(allInfoLength);
	vector<int> codeBits(allBlockLength);


	for (int d = 0; d < DEVISION; d++)
	{
		vector<int> subInfoBits(infoLength);
		for (int b = 0; b < infoLength; b++)
			subInfoBits.at(b) = rand() % 2;


		copy(subInfoBits.begin(), subInfoBits.end(), infoBits.begin() + (infoLength * d));

		vector<int> subCodeBits = SubBlockEncoding(subInfoBits);

		copy(subCodeBits.begin(), subCodeBits.end(), codeBits.begin() + (blockLength * d));
	}


	return codeBits;
}

vector<int> LDPC::MutlipleSubBlockDecodeer(vector<double> receivedLLR, vector<double> &decodedLLR, int &totalErr , int times)
{
	vector<double> LLR(receivedLLR.size());
	copy(receivedLLR.begin(), receivedLLR.end(), LLR.begin());


	vector<int> decodedCodeword(receivedLLR.size());	decodedLLR.resize(receivedLLR.size());
	vector<bool> preCheckTable(DEVISION), sufCheckTable(DEVISION);	vector<int> preErrTable(DEVISION), sufErrTable(DEVISION);
	
	if (times == 0) 
	{
		for (int d = 0; d < DEVISION; d++)
		{
			preCheckTable.at(d) = 0;
			preErrTable.at(d) = 0;
		}
	}
	else 
	{
		copy(checkTable[times - 1].begin(), checkTable[times - 1].end(), preCheckTable.begin());
		copy(errTable[times - 1].begin(), errTable[times - 1].end(), preErrTable.begin());
	}

	for (int d = 0; d < DEVISION; d++)
	{
		vector<double> subLLR(blockLength);
		copy(LLR.begin() + (d*blockLength), LLR.begin() + ((d + 1)*blockLength), subLLR.begin());
		if (preCheckTable.at(d) != true)
		{
			bool check = 0;	int errNum = 0;
			vector<double> subDecodedLLR;
			vector<int> subCodeword = SubBlockDecoding(subLLR, subDecodedLLR, d, errNum, check);
			copy(subCodeword.begin(), subCodeword.end(), decodedCodeword.begin() + (d*blockLength));
			copy(subDecodedLLR.begin(), subDecodedLLR.end(), decodedLLR.begin() + (d*blockLength));
			sufCheckTable.at(d) = check;	sufErrTable.at(d) = errNum;
		}
		else 
		{
			vector<double> subDecodedLLR(subLLR.size());
			vector<int> subCodeword(subLLR.size());

			for (int j = 0; j < subLLR.size(); j++)
			{
				subDecodedLLR.at(j) = subLLR.at(j);
				(subLLR.at(j) >= 0) ? subCodeword.at(j) = 0 : subCodeword.at(j) = 1;
			}

			copy(subCodeword.begin(), subCodeword.end(), decodedCodeword.begin() + (d*blockLength));
			copy(subDecodedLLR.begin(), subDecodedLLR.end(), decodedLLR.begin() + (d*blockLength));

			sufCheckTable.at(d) = preCheckTable.at(d);	sufErrTable.at(d) = preErrTable.at(d);
		}
		
	}

	int sumCheck = 0, sumErr = 0;
	for (int j = 0; j < sufCheckTable.size(); j++)
	{
		sumCheck += sufCheckTable.at(j);
		sumErr += sufErrTable.at(j);
	}

	checkTable.push_back(sufCheckTable);
	errTable.push_back(sufErrTable);
	errFlag.push_back (sumCheck);
	totalErr = sumErr;

	return decodedCodeword;
}


