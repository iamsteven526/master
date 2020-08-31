#ifndef MATH_OPERATION
#define MATH_OPERATION
#define _USE_MATH_DEFINES
#include"setting.h"
#include<math.h>
#include <cstring>
#include<malloc.h>

class FFT
{
public:
	FFT(int _size, int _layer) : FFT_size(_size), layer(_layer)
	{
		bitReversalOrder.resize(FFT_size);
		create_bit_rev_order();
		j = complexNum(0, 1);
		omega = exp(complexNum(-2 * M_PI / FFT_size, 0)*j);

	}

	vector<complexNum> FastFourierTransform(vector<complexNum> inputSequence)
	{
		if (inputSequence.size() % FFT_size != 0)
		{
			int lack = FFT_size - (inputSequence.size() % FFT_size);
			for (int i = 0; i < lack; i++)
				inputSequence.push_back(complexNum(0, 0));
		}

		vector <complexNum> outputSequence(inputSequence.size());
		vector<complexNum> bitReversalFFTinput(FFT_size);

		for (int k = 0; k < FFT_size; k++)
			bitReversalFFTinput.at(k) = inputSequence.at(bitReversalOrder.at(k));

		vector<complexNum> tmp(FFT_size);


		for (int l = 0; l < layer; l++)
		{
			int subsize = int(pow(2, l + 1));
			int count = FFT_size / subsize;

			for (int c = 0; c < count; c++)
			{
				for (int s = 0; s < subsize / 2; s++)
				{
					tmp.at(c*subsize + s) = bitReversalFFTinput.at(c*subsize + s) + bitReversalFFTinput.at(c*subsize + s + pow(2, l))*pow(omega, s*count);
					tmp.at(c*subsize + s + pow(2, l)) = bitReversalFFTinput.at(c*subsize + s) - bitReversalFFTinput.at(c*subsize + s + pow(2, l))*pow(omega, s*count);
				}
			}
			copy(tmp.begin(), tmp.end(), bitReversalFFTinput.begin());
		}

		copy(bitReversalFFTinput.begin(), bitReversalFFTinput.end(), outputSequence.begin());


		return outputSequence;
	}

	vector<complexNum> InverseFastFourierTransform(vector<complexNum> inputSequence)
	{
		if (inputSequence.size() % FFT_size != 0)
		{
			for (int i = 0; i < FFT_size - (inputSequence.size() % FFT_size); i++)
				inputSequence.push_back(complexNum(0, 0));
		}

		vector <complexNum> outputSequence(inputSequence.size());
		vector<complexNum> conjIFFTinput(FFT_size);
		vector<complexNum> FFTconjIFFTinput(FFT_size);


		for (int k = 0; k < FFT_size; k++)
			conjIFFTinput.at(k) = conj(inputSequence.at(k));

		FFTconjIFFTinput = FastFourierTransform(conjIFFTinput);

		for (int k = 0; k < FFT_size; k++)
			outputSequence.at(k) = conj(FFTconjIFFTinput.at(k)) / complexNum(FFT_size, 0);


		return outputSequence;
	}

	vector<complexNum> NormalizedFastFourierTransform(vector<complexNum> inputSequence)				//normalized
	{
		if (inputSequence.size() % FFT_size != 0)
		{
			int lack = FFT_size - (inputSequence.size() % FFT_size);
			for (int i = 0; i < lack; i++)
				inputSequence.push_back(complexNum(0, 0));
		}

		vector <complexNum> outputSequence(inputSequence.size());
		vector<complexNum> bitReversalFFTinput(FFT_size);

		for (int k = 0; k < FFT_size; k++)
			bitReversalFFTinput.at(k) = inputSequence.at(bitReversalOrder.at(k));

		vector<complexNum> tmp(FFT_size);


		for (int l = 0; l < layer; l++)
		{
			int subsize = int(pow(2, l + 1));
			int count = FFT_size / subsize;

			for (int c = 0; c < count; c++)
			{
				for (int s = 0; s < subsize / 2; s++)
				{
					tmp.at(c*subsize + s) = bitReversalFFTinput.at(c*subsize + s) + bitReversalFFTinput.at(c*subsize + s + pow(2, l))*pow(omega, s*count);
					tmp.at(c*subsize + s + pow(2, l)) = bitReversalFFTinput.at(c*subsize + s) - bitReversalFFTinput.at(c*subsize + s + pow(2, l))*pow(omega, s*count);
				}
			}
			copy(tmp.begin(), tmp.end(), bitReversalFFTinput.begin());
		}

		copy(bitReversalFFTinput.begin(), bitReversalFFTinput.end(), outputSequence.begin());


		for (int k = 0; k < outputSequence.size(); k++)
			outputSequence.at(k) = complexNum(1. / sqrt(FFT_size), 0)*outputSequence.at(k);

		return outputSequence;
	}

	vector<complexNum> NormalizedInverseFastFourierTransform(vector<complexNum> inputSequence)		//normalized
	{
		if (inputSequence.size() % FFT_size != 0)
		{
			for (int i = 0; i < FFT_size - (inputSequence.size() % FFT_size); i++)
				inputSequence.push_back(complexNum(0, 0));
		}

		vector <complexNum> outputSequence(inputSequence.size());
		vector<complexNum> conjIFFTinput(FFT_size);
		vector<complexNum> FFTconjIFFTinput(FFT_size);


		for (int k = 0; k < FFT_size; k++)
			conjIFFTinput.at(k) = conj(inputSequence.at(k));

		FFTconjIFFTinput = FastFourierTransform(conjIFFTinput);

		for (int k = 0; k < FFT_size; k++)
			outputSequence.at(k) = conj(FFTconjIFFTinput.at(k)) / complexNum(FFT_size, 0);

		for (int k = 0; k < outputSequence.size(); k++)
			outputSequence.at(k) = complexNum(sqrt(FFT_size), 0)*outputSequence.at(k);

		return outputSequence;
	}

private:

	complexNum j;
	complexNum omega;
	int FFT_size, layer;

	vector<uint16_t> bitReversalOrder;

	void create_bit_rev_order() {																				//bit reversal
		for (uint16_t i = 0; i < FFT_size; ++i) {
			uint16_t to_be_reversed = i;
			bitReversalOrder.at(i) = (uint16_t)((to_be_reversed & 1) << (layer - 1));
			for (uint8_t j = (uint8_t)(layer - 1); j; --j) {
				to_be_reversed >>= 1;
				bitReversalOrder.at(i) += (to_be_reversed & 1) << (j - 1);
			}
		}
	}
};

class COMPLEX_OPERATION
{
public:
	complexNum** CreateMatrix(int row, int col)
	{
		complexNum **p = new complexNum *[row];
		for (int i = 0; i < row; i++)
		{
			p[i] = new complexNum[col];
			memset(p[i], 0, col * sizeof(complexNum));
		}
		return p;
	};

	void DeleteMatrix(complexNum **matrix)
	{
		int row = malloc_usable_size(matrix) / sizeof(matrix[0]);
		for (int i = 0; i < row; i++)
			delete[] matrix[i];
		delete[] matrix;
	}

	complexNum** ComplexInverseMatrix(complexNum **matrix)
	{
		int row = malloc_usable_size(matrix) / sizeof(matrix[0]);
		int col = malloc_usable_size(matrix[0]) / sizeof(matrix[0][0]);

		if (row != col)
			throw length_error("#error: This matrix is not a square matrix. ");

		complexNum **tmpArray = CreateMatrix(row, 2 * col);

		for (int i = 0; i < row; i++)
		{
			tmpArray[i][i + col] = complexNum(1, 0);
			for (int j = 0; j < col; j++)
				tmpArray[i][j] = matrix[i][j];
		}

		for (int r = 0; r < row; r++)
		{
			for (int idx = r; idx < row; idx++)
			{
				if (tmpArray[idx][r] != complexNum(0, 0))
				{
					for (int col_idx = 0; col_idx < 2 * col; col_idx++)
						swap(tmpArray[r][col_idx], tmpArray[idx][col_idx]);
					break;
				}
			}

			complexNum scalar = tmpArray[r][r];
			for (int c = 0; c < 2 * col; c++)
				tmpArray[r][c] = tmpArray[r][c] / scalar;

			for (int idx = 0; idx < row; idx++)
			{
				if (idx != r && tmpArray[idx][r] != complexNum(0, 0))
				{
					complexNum Mul = tmpArray[idx][r];
					for (int c = 0; c < 2 * col; c++)
					{
						tmpArray[idx][c] = tmpArray[idx][c] - Mul * tmpArray[r][c];
					}
				}
			}

		}

		complexNum **invMatrix = CreateMatrix(row, col);
		for (int i = 0; i < row; i++)
		{
			for (int j = 0; j < col; j++)
				invMatrix[i][j] = tmpArray[i][j + col];
		}

		DeleteMatrix(tmpArray);

		return invMatrix;

	}

	inline complexNum** HermitianOperator(complexNum **matrix)
	{
		int row = malloc_usable_size(matrix) / sizeof(matrix[0]);
		int col = malloc_usable_size(matrix[0]) / sizeof(matrix[0][0]);
		complexNum** matrixHermitian = CreateMatrix(col, row);

		for (int r = 0; r < row; r++)
		{
			for (int c = 0; c < col; c++)
				matrixHermitian[c][r] = complexNum(matrix[r][c].real(), (-1)* matrix[r][c].imag());
		}

		return matrixHermitian;
	}

	complexNum** MatrixMultiplication(complexNum **leftMatrix, complexNum **rightMatrix)
	{
		int rowLeft = malloc_usable_size(leftMatrix) / sizeof(leftMatrix[0]);
		int colLeft = malloc_usable_size(leftMatrix[0]) / sizeof(leftMatrix[0][0]);

		int rowRight = malloc_usable_size(rightMatrix) / sizeof(rightMatrix[0]);
		int colRight = malloc_usable_size(rightMatrix[0]) / sizeof(rightMatrix[0][0]);

		if (colLeft != rowRight)
			throw length_error("#error: Column length of the left matrix isn't equal to the row length of the right matrix ");

		complexNum **matrix = CreateMatrix(rowLeft, colRight);

		for (int i = 0; i < rowLeft; i++)
		{
			for (int j = 0; j < colRight; j++)
			{
				for (int idx = 0; idx < rowRight; idx++)
					matrix[i][j] += (leftMatrix[i][idx] * rightMatrix[idx][j]);
			}
		}

		return matrix;
	}

private:
	double** RealCreateMatrix(int row, int col)
	{
		double **p = new double *[row];
		for (int i = 0; i < row; i++)
		{
			p[i] = new double[col];
			memset(p[i], 0, col * sizeof(double));
		}
		return p;
	};

	void RealDeleteMatrix(double **matrix)
	{
		int row = malloc_usable_size(matrix) / sizeof(matrix[0]);
		for (int i = 0; i < row; i++)
			delete[] matrix[i];
		delete[] matrix;
	}

	double** RealMatrixMultiplication(double **leftMatrix, double **rightMatrix)
	{
		int rowLeft = malloc_usable_size(leftMatrix) / sizeof(leftMatrix[0]);
		int colLeft = malloc_usable_size(leftMatrix[0]) / sizeof(leftMatrix[0][0]);

		int rowRight = malloc_usable_size(rightMatrix) / sizeof(rightMatrix[0]);
		int colRight = malloc_usable_size(rightMatrix[0]) / sizeof(rightMatrix[0][0]);

		if (colLeft != rowRight)
			throw length_error("#error: Column length of the left matrix isn't equal to the row length of the right matrix ");

		double **matrix = RealCreateMatrix(rowLeft, colRight);

		for (int i = 0; i < rowLeft; i++)
		{
			for (int j = 0; j < colRight; j++)
			{
				for (int idx = 0; idx < rowRight; idx++)
					matrix[i][j] += (leftMatrix[i][idx] * rightMatrix[idx][j]);
			}
		}

		return matrix;
	}

	double** RealInverseMatix(double **matrix)
	{
		int row = malloc_usable_size(matrix) / sizeof(matrix[0]);
		double detMatrix = Determination(matrix);
		if (detMatrix == 0)
		{
			cout << "not invertible!" << endl;
			throw length_error("#error: not invertible! ");
		}

		double** coefficentMatrix = RealCreateMatrix(row, row);
		double** realIverseMatrix = RealCreateMatrix(row, row);

		for (int r = 0; r < row; r++)
		{
			for (int c = 0; c < row; c++)
			{
				double **subMatrix = RealCreateMatrix(row - 1, row - 1);

				int subRow = 0;
				for (int m = 0; m < row; m++)
				{
					if (m != r)
					{
						int subCol = 0;
						for (int n = 0; n < row; n++)
						{
							if (n != c)
							{
								subMatrix[subRow][subCol] = matrix[m][n];
								subCol++;
							}
						}
						subRow++;
					}
				}


				coefficentMatrix[r][c] = (1. / detMatrix)*pow(-1, (r + 1) + (c + 1))*Determination(subMatrix);

				RealDeleteMatrix(subMatrix);
			}
		}

		for (int r = 0; r < row; r++)
		{
			for (int c = 0; c < row; c++)
			{
				realIverseMatrix[c][r] = coefficentMatrix[r][c];
				if (isnan(norm(realIverseMatrix[r][c])) == 1)
				{
					cout << endl << "det: " << detMatrix << endl;
					cout << endl << "matrix: " << endl;
					for (int m = 0; m < row; m++)
					{
						for (int n = 0; n < row; n++)
							cout << matrix[m][n] << " ";
						cout << endl;
					}
					cout << endl << "invs matrix: " << endl;
					for (int m = 0; m < row; m++)
					{
						for (int n = 0; n < row; n++)
							cout << realIverseMatrix[m][n] << " ";
						cout << endl;
					}
					break; break;
				}
			}
		}
		RealDeleteMatrix(coefficentMatrix);

		return realIverseMatrix;
	}

	double Determination(double **matrix)				// Barriess Algorithm
	{
		int row = malloc_usable_size(matrix) / sizeof(matrix[0]);

		double **matrix1 = RealCreateMatrix(row, row);
		for (int r = 0; r < row; r++)
			memcpy(matrix1[r], matrix[r], row * sizeof(double));

		if (matrix1[0][0] == 0)
		{
			for (int r = 1; r < row; r++)
			{
				if (matrix1[r][0] != 0)
				{
					for (int c = 0; c < row; c++)
					{
						matrix1[0][c] = (-1)* matrix1[0][c];
						swap(matrix1[0][c], matrix1[r][c]);
					}
					break;
				}
			}
		}


		double p = 1;
		double **tmpMatrix = RealCreateMatrix(row, row);
		for (int idx = 0; idx < row - 1; idx++)
		{
			for (int r = 0; r < row; r++)
				memcpy(tmpMatrix[r], matrix1[r], row * sizeof(double));

			for (int r = idx + 1; r < row; r++)
			{
				for (int c = 0; c < row; c++)
					matrix1[r][c] = (tmpMatrix[idx][idx] * tmpMatrix[r][c] - tmpMatrix[r][idx] * tmpMatrix[idx][c]) / p;
			}
			p = matrix1[idx][idx];

		}
		double determineValue = matrix1[row - 1][row - 1];

		RealDeleteMatrix(matrix1);
		RealDeleteMatrix(tmpMatrix);

		return determineValue;
	}

};
#endif // !MATH_OPERATION

