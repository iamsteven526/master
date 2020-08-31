#include"channel_coding.h"
#include<stdexcept>


#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) > (y) ? (x) : (y))

int** BitMatrix::CreateMatrix(int row, int col)
{
	int **p = new int *[row];
	for (int i = 0; i < row; i++)
	{
		p[i] = new int[col];
		memset(p[i], 0, col * sizeof(int));
	}
	return p;
}

void  BitMatrix::DeleteMatrix(int **matrix)
{
	int row = _msize(matrix) / sizeof(matrix[0]);
	for (int i = 0; i < row; i++)
		delete[] matrix[i];
	delete[] matrix;
}

int** BitMatrix::MatrixMultiplication(int **matrix_left, int **matrix_right)
{
	int row_left = _msize(matrix_left) / sizeof(matrix_left[0]);
	int col_left = _msize(matrix_left[0]) / sizeof(matrix_left[0][0]);

	int row_right = _msize(matrix_right) / sizeof(matrix_right[0]);
	int col_right = _msize(matrix_right[0]) / sizeof(matrix_right[0][0]);

	if (col_left != row_right)
		throw length_error("#error: Column length of the left matrix isn't equal to the row length of the right matrix ");
	
	int **matrix = CreateMatrix(row_left, col_right);

	for (int i = 0; i < row_left; i++) 
	{
		for (int j = 0; j < col_right; j++) 
		{
			for (int idx = 0; idx < row_right; idx++)
				matrix[i][j] ^= (matrix_left[i][idx] & matrix_right[idx][j]);
		}
	}
	return matrix;
}

int BitMatrix::RankCalculation(int **matrix)
{
	int row = _msize(matrix) / sizeof(matrix[0]); 
	int col = _msize(matrix[0]) / sizeof(matrix[0][0]); 

	int length = min(row, col);

	int **tmp_matrix = CreateMatrix(row, col);

	for (int i = 0; i < row; i++)
		memcpy(tmp_matrix[i], matrix[i], col * sizeof(int));

	for (int row_idx = 0; row_idx < length; row_idx++)
	{
		for (int idx = row_idx; idx < row; idx++)
		{
			if (tmp_matrix[idx][row_idx] == 1)
			{
				for (int col_idx = 0; col_idx < col; col_idx++)
					swap(tmp_matrix[row_idx][col_idx], tmp_matrix[idx][col_idx]);
				break;
			}
		}
		for (int idx = 0; idx < row; idx++)
		{
			if (idx != row_idx && tmp_matrix[idx][row_idx] == 1) {
				for (int col_idx = 0; col_idx < col; col_idx++)
				{
					tmp_matrix[idx][col_idx] = tmp_matrix[idx][col_idx] ^ tmp_matrix[row_idx][col_idx];
				}
			}
		}
	}

	int rank = row;

	for (int i = 0; i < row; i++) 
	{
		int sum = 0;
		for (int j = 0; j < col; j++) 
			sum += tmp_matrix[i][j];
		if (sum == 0)
			rank--;
	}
	DeleteMatrix(tmp_matrix);
	return rank;
}

int** BitMatrix::InverseMatrix(int **matrix)
{
	int row = _msize(matrix) / sizeof(matrix[0]);
	int col = _msize(matrix[0]) / sizeof(matrix[0][0]);

	if (row != col)
		throw length_error("#error: This matrix is not a square matrix. ");

	int **tmp_arr = CreateMatrix(row, 2 * col);
	for (int i = 0; i < row; i++)
	{
		tmp_arr[i][i + col] = 1;
		for (int j = 0; j < col; j++)
			tmp_arr[i][j ] = matrix[i][j];
	}


	for (int row_idx = 0; row_idx < row; row_idx++)
	{
		for (int idx = row_idx; idx < row; idx++)
		{
			if (tmp_arr[idx][row_idx] == 1)
			{
				for (int col_idx = 0; col_idx < 2 * col; col_idx++)
					swap(tmp_arr[row_idx][col_idx], tmp_arr[idx][col_idx]);
				break;
			}
		}
		for (int idx = 0; idx < row; idx++)
		{
			if (idx != row_idx && tmp_arr[idx][row_idx] == 1) {
				for (int col_idx = 0; col_idx < 2*col; col_idx++)
				{
					tmp_arr[idx][col_idx] = tmp_arr[idx][col_idx] ^ tmp_arr[row_idx][col_idx];
				}
			}
		}
	}
	bool check_identity = true;
	for (int i = 0; i < row; i++) 
	{
		if (tmp_arr[i][i] != 1) 
		{
			check_identity = false;
			break;
		}
		for (int j = 0; j < col; j++){
			if (i != j){
				if (tmp_arr[i][j] != 0){
					check_identity = false;
					break;
				}
			}
		}
	}

	//cout << "Inversible check: " << check_identity << endl;

	int **inv_matrix = CreateMatrix(row, col);
	for (int i = 0; i < row; i++) 
	{
		for (int j = 0; j < col; j++)
			inv_matrix[i][j] = tmp_arr[i][j + col];
	}

	DeleteMatrix(tmp_arr);

	return inv_matrix;

}

vector<uint16_t> BitMatrix::FindIndependentRow(int **matrix)
{
	int row = _msize(matrix) / sizeof(matrix[0]);
	int col = _msize(matrix[0]) / sizeof(matrix[0][0]);
	int length = min(row, col);

	vector<uint16_t> row_operation_order(row);
	for (int i = 0; i < row; i++)
		row_operation_order.at(i) = i;

	int **tmp_matrix = CreateMatrix(row, col);

	for (int i = 0; i < row; i++)
		memcpy(tmp_matrix[i], matrix[i], col * sizeof(int));

	for (int row_idx = 0; row_idx < length; row_idx++)
	{
		for (int idx = row_idx; idx < row; idx++)
		{
			if (tmp_matrix[idx][row_idx] == 1)
			{
				swap(row_operation_order.at(row_idx), row_operation_order.at(idx));
				for (int col_idx = 0; col_idx < col; col_idx++)
					swap(tmp_matrix[row_idx][col_idx], tmp_matrix[idx][col_idx]);
				break;
			}
		}

		for (int idx = 0; idx < row; idx++)
		{
			if (idx != row_idx && tmp_matrix[idx][row_idx] == 1) {
				for (int col_idx = 0; col_idx < col; col_idx++)
					tmp_matrix[idx][col_idx] = tmp_matrix[idx][col_idx] ^ tmp_matrix[row_idx][col_idx];
			}
		}

	}

	
	DeleteMatrix(tmp_matrix);
	return row_operation_order;

}
