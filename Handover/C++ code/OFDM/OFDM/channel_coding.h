#ifndef CHANNEL_CODING
#define CHANNEL_CODING
#include"setting.h"
#include <cstring>

class BitMatrix 
{
public:
	int** CreateMatrix(int row, int col);
	void DeleteMatrix(int **matrix);
	int** MatrixMultiplication(int **matrix_left, int **matrix_right);
	int RankCalculation(int **matrix);
	int** InverseMatrix(int **matrix);
	vector<uint16_t> FindIndependentRow(int **matrix);

};

class LDPC 
{
public:
	LDPC()
	{
		blockLength = MOD_BIT * (N / DEVISION);
		infoLength = blockLength * R;
		allBlockLength = MOD_BIT * N;
		allInfoLength = allBlockLength * R;

		if (R != 1) 
		{
			H = new int *[blockLength - infoLength];
			for (int i = 0; i < blockLength - infoLength; i++)
			{
				H[i] = new int[blockLength];
				memset(H[i], 0, blockLength * sizeof(int));
			}


			FILE *Hfile;
			Hfile = fopen(H_FILE, "r");
			for (int row_idx = 0; row_idx < blockLength - infoLength; row_idx++)
			{
				for (int col_idx = 0; col_idx < blockLength; col_idx++)
					fscanf(Hfile, "%d", &H[row_idx][col_idx]);
			}fclose(Hfile);

			checkNodeSet.resize(0, vector<int>(0));
			for (int col = 0; col < blockLength; col++)
			{
				vector<int> oneVariable2CheckSet;
				for (int row = 0; row < (blockLength - infoLength); row++)
				{
					if (H[row][col] == 1)
						oneVariable2CheckSet.push_back(row);
				}
				checkNodeSet.push_back(oneVariable2CheckSet);
			}

			variableNodeSet.resize(0, vector<int>(0));
			for (int row = 0; row < (blockLength - infoLength); row++)
			{
				vector<int> oneCheck2VariableSet;
				for (int col = 0; col < blockLength; col++)
				{
					if (H[row][col] == 1)
						oneCheck2VariableSet.push_back(col);
				}
				variableNodeSet.push_back(oneCheck2VariableSet);
			}

			G_P = CreateGeneratorMatrix(rowOperationOrder);
			colOperationOrder = CreateRow2colOrder(rowOperationOrder);
		}
	}

	int **CreateGeneratorMatrix(vector<uint16_t> &row_operation_order);

	vector<int> SubBlockEncoding(vector<int> subInfoBits);
	vector<int> SubBlockDecoding(vector<double> receivedLLR, vector<double> &subDecodedLLR,int index, int &err, bool &check);

	vector<uint16_t> CreateRow2colOrder(vector<uint16_t> row_operation_order);
	vector<int> MutlipleSubBlockEncoder();
	vector<int> MutlipleSubBlockDecodeer(vector<double> receivedLLR, vector<double> &decodedLLR,int &totalErr, int times);

	vector<vector<bool>> checkTable;
	vector<vector<int>> errTable;
	vector<int> errFlag;

	vector<vector<int>> LSdecodeCodeBits;
	vector<vector<double>> LSdecodedLLR;

private:
	int **H, **G_P;
	vector<int> infoBits;
	
	int blockLength, infoLength, allBlockLength, allInfoLength;
	vector<uint16_t> colOperationOrder, rowOperationOrder;

	vector<vector<int>> checkNodeSet;   // for variable node calculation
	vector<vector<int>> variableNodeSet; // for check node calculation

};
#endif
