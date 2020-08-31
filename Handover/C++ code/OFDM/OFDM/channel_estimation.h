#ifndef CHANNEL_ESTIMATION
#define CHANNEL_ESTIMATION
#include"setting.h"
#include"math_operation.h"
#include"modulation.h"

COMPLEX_OPERATION complexOperation;
QAM qam(MOD_BIT);
FFT fft(N, LAYER);

class ESTIMATION
{
public:

	vector<complexNum> Permutation(vector<complexNum> X)
	{
		vector<complexNum> output(X.size());

		vector<vector<complexNum>> Xarray(N / S, vector<complexNum>(S));

		for (int r = 0; r < N / S; r++)
			copy(X.begin() + (r*S), X.begin() + ((r + 1)*S), Xarray[r].begin());
		
		segmentSwapIdx.resize(N / S);
		int sIdx = 0;
		for (int s = 0; s < (N / (S*DEVISION)); s++)
		{
			for (int d = 0; d < DEVISION; d++)
			{
				copy(Xarray[d*(N / (S*DEVISION)) + s].begin(), Xarray[d*(N / (S*DEVISION)) + s].end(), output.begin() + (sIdx * S));
				segmentSwapIdx.at(sIdx) = d * (N / (S*DEVISION)) + s;
				sIdx++;
				
			}
		}

		return output;
	}

	vector<complexNum> InversePermutation(vector<complexNum> Y, vector<double> CFR, vector<double> &permutatedCFR)
	{
		vector<complexNum> output(Y.size());
		permutatedCFR.resize(CFR.size());
	
		vector<vector<complexNum>> Yarray(N / S, vector<complexNum>(S));
		vector<vector<double>> CFRarray(N / S, vector<double>(S));
		for (int r = 0; r < N / S; r++)
		{
			copy(Y.begin() + (r*S), Y.begin() + ((r + 1)*S), Yarray[r].begin());
			copy(CFR.begin() + (r*S), CFR.begin() + ((r + 1)*S), CFRarray[r].begin());
		}
		
		int sIdx = 0;
		for (int d = 0; d < DEVISION; d++) 
		{
			for (int s = 0; s < (N / (S*DEVISION)); s++) 
			{
				copy(Yarray[sIdx].begin(), Yarray[sIdx].end(), output.begin() + (segmentSwapIdx.at(sIdx)*S));
				copy(CFRarray[sIdx].begin(), CFRarray[sIdx].end(), permutatedCFR.begin() + (segmentSwapIdx.at(sIdx)*S));
				sIdx++;
			}
		}

		return output;
	}

	complexNum** CreateHmatrix(vector<vector<complexNum>> timeVaryingCIR, int startIdx)
	{
		complexNum **H = complexOperation.CreateMatrix(N, N);

		complexNum **channelMatrix = complexOperation.CreateMatrix(N, N);
		complexNum **leftMatrix = complexOperation.CreateMatrix(N, N);

		for (int time = 0; time < N; time++)
		{
			for (int l = 0; l < L; l++)
				channelMatrix[time][(time - l + N) % N] = timeVaryingCIR.at(l).at(time + startIdx);
		}

		for (int c = 0; c < N; c++)
		{
			vector<complexNum> input(N);
			for (int i = 0; i < N; i++)
				input.at(i) = channelMatrix[i][c];

			vector<complexNum> subColumn = fft.NormalizedFastFourierTransform(input);
			for (int i = 0; i < N; i++)
				leftMatrix[i][c] = subColumn.at(i);
		}

		for (int r = 0; r < N; r++)
		{
			vector<complexNum> input(N);
			for (int i = 0; i < N; i++)
				input.at(i) = leftMatrix[r][i];

			vector<complexNum> subColumn = fft.NormalizedInverseFastFourierTransform(input);
			for (int i = 0; i < N; i++)
				H[r][i] = subColumn.at(i);

		}

		complexOperation.DeleteMatrix(channelMatrix);
		complexOperation.DeleteMatrix(leftMatrix);

		return H;
	}

	vector<complexNum>  OneTapEqualizer(vector<complexNum> CFR, vector<complexNum> Y)
	{
		vector<complexNum> estimatedOutput(N);
		for (int i = 0; i < N; i++)
			estimatedOutput.at(i) = Y.at(i) / CFR.at(i);

		return estimatedOutput;
	}

	vector<complexNum> PartialMMSEequalizer(complexNum **H, vector<complexNum> Y, double N0 ,vector<double> &channelPower)
	{
		complexNum **sparseH = complexOperation.CreateMatrix(N, N);
		for (int d = 0; d < N; d++)
		{
			for (int k = 0; k < N; k++)
			{
				if (abs(d - k) <= D || abs(d - k) >= (N - D))
					sparseH[d][k] = H[d][k];
			}
		}

		channelPower.resize(N);
		vector<complexNum> estimated_X(N, 0);
		for (int k = 0; k < N; k++)
		{
			complexNum **Hk = complexOperation.CreateMatrix(2 * D + 1, 4 * D + 1);

			for (int r = 0; r <= 2 * D; r++)
			{
				for (int c = 0; c <= 4 * D; c++)
				{
					int row_idx = (k - D + r + N) % N;
					int col_idx = (k - 2 * D + c + N) % N;
					Hk[r][c] = sparseH[row_idx][col_idx]; // +N for mod
				}
			}

			vector<complexNum> Yk(2 * D + 1);
			for (int i = 0; i < 2 * D + 1; i++)
				Yk.at(i) = Y.at((k - D + i + N) % N);			//+N for mod

			estimated_X.at(k) = PartialMMSEequalization(Hk, Yk, N0, channelPower.at(k));

			complexOperation.DeleteMatrix(Hk);
			
		}

		complexOperation.DeleteMatrix(sparseH);

		return estimated_X;
	}

	vector<complexNum> LMMSEequalizer(complexNum **H, vector<complexNum> Y, double N0 ,vector<double> &channelPower) 
	{
		complexNum **Her_H = complexOperation.HermitianOperator(H);
		complexNum **Mul_H_Her_H = complexOperation.MatrixMultiplication(H, Her_H);
		for (int i = 0; i < N; i++)
			Mul_H_Her_H[i][i] += complexNum(N0, 0);

		complexNum **inv_mul_H_Her_H = complexOperation.ComplexInverseMatrix(Mul_H_Her_H);
		complexNum **Her_G = complexOperation.MatrixMultiplication(Her_H, inv_mul_H_Her_H);

		complexOperation.DeleteMatrix(Her_H);
		complexOperation.DeleteMatrix(Mul_H_Her_H);
		complexOperation.DeleteMatrix(inv_mul_H_Her_H);

		vector<complexNum> estimated_X(N, 0);
		for (int i = 0; i < estimated_X.size(); i++)
		{
			for (int j = 0; j < Y.size(); j++)
				estimated_X.at(i) += Her_G[i][j] * Y.at(j);
		}

		channelPower.resize(N);
		
		for (int j = 0; j < N; j++)
		{
			double sum = 0;

			for (int r = 0; r < N; r++)
				sum += abs(Her_G[j][r]*conj(Her_G[j][r]));

			channelPower.at(j) = sum;

		}

		complexOperation.DeleteMatrix(Her_G);

		return estimated_X;
	}

	vector<complexNum> PartialMMSESIC(complexNum **H, vector<complexNum> Y, double N0)
	{
		
		complexNum **sparseH = complexOperation.CreateMatrix(N, N);
		for (int r = 0; r < N; r++)
		{
			for (int c = 0; c < N; c++)
			{
				if (abs(r - c) <= D || abs(r - c) >= (N - D))
					sparseH[r][c] = H[r][c];
			}
		}

		vector<complexNum> Ycopy(Y.size());
		copy(Y.begin(), Y.end(), Ycopy.begin());

		int check = 1;
		vector<bool> state(N, true);
		vector<complexNum> estimated_S(N);
		vector<complexNum> estimated_X(N, 0);

		vector<double> diagonalPower(N, 0);

		while (check != 0)
		{
			for (int d = 0; d < N; d++)
				diagonalPower.at(d) = norm(sparseH[d][d]);

			int k = max_element(diagonalPower.begin(), diagonalPower.end()) - diagonalPower.begin();

			complexNum **Hk = complexOperation.CreateMatrix(2 * D + 1, 4 * D + 1);
			for (int r = 0; r < 2 * D + 1; r++)
			{
				for (int c = 0; c < 4 * D + 1; c++)
				{
					int row_idx = (k - D + r + N) % N;
					int col_idx = (k - 2 * D + c + N) % N;
					Hk[r][c] = sparseH[row_idx][col_idx]; // +N for mod
				}
			}

			vector<complexNum> Yk(2 * D + 1);
			for (int i = 0; i < 2 * D + 1; i++)
				Yk.at(i) = Ycopy.at((k - D + i + N) % N);		//+N for mod

			double gainPower = 0;
			estimated_X.at(k) = PartialMMSEequalization(Hk, Yk, N0, gainPower);

			complexOperation.DeleteMatrix(Hk);
			
			estimated_S.at(k) = qam.oneSampleHardDecision(estimated_X.at(k));

			for (int j = 0; j < N; j++)
			{
				Ycopy.at(j) = Ycopy.at(j) - sparseH[j][k] * estimated_S.at(k);
				sparseH[j][k] = complexNum(0, 0);
			}

			state.at(k) = 0;

			int checkSum = 0;
			for (int r = 0; r < N; r++)
				checkSum += state.at(r);
			check = checkSum;
		}

		complexOperation.DeleteMatrix(sparseH);

		return estimated_X;
	}

	vector<complexNum> ListPartialMMSESIC(complexNum **H, vector<complexNum> Y, double N0)
	{
		vector<complexNum> Ycopy(Y.size());
		copy(Y.begin(), Y.end(), Ycopy.begin());
		vector<complexNum> output(N);

		vector<vector<complexNum>> path, Youtput;
		path.push_back(vector<complexNum>(N, 0));
		Youtput.push_back(vector<complexNum>(N, 0));
		for (int c = 0; c < N; c++)
			Youtput[0].at(c) = Ycopy.at(c);

		int check = 1;
		vector<bool> state(N, true);
		vector<double> diagonalPower(N, 0), probablity(1, 1);

		complexNum **sparseH = complexOperation.CreateMatrix(N, N);
		for (int r = 0; r < N; r++)
		{
			for (int c = 0; c < N; c++)
			{
				if (abs(r - c) <= D || abs(r - c) >= (N - D))
					sparseH[r][c] = H[r][c];
			}
		}

		while (check != 0)
		{
			for (int d = 0; d < N; d++)
				diagonalPower.at(d) = norm(sparseH[d][d]);

			int k = max_element(diagonalPower.begin(), diagonalPower.end()) - diagonalPower.begin();

			complexNum **Hk = complexOperation.CreateMatrix(2 * D + 1, 4 * D + 1);
			for (int r = 0; r < 2 * D + 1; r++)
			{
				for (int c = 0; c < 4 * D + 1; c++)
				{
					int row_idx = (k - D + r + N) % N;
					int col_idx = (k - 2 * D + c + N) % N;
					Hk[r][c] = sparseH[row_idx][col_idx]; // +N for mod
				}
			}

			vector<vector<complexNum>> bufferPath, bufferY;
			vector<double> bufferProb;

			vector<complexNum> Yk(2 * D + 1);

			for (int list = 0; list < path.size(); list++)
			{
				for (int i = 0; i < 2 * D + 1; i++)
					Yk.at(i) = Youtput[list].at((k - D + i + N) % N);

				double gainPower;
				complexNum valuek = PartialMMSEequalization(Hk, Yk, N0, gainPower);

				vector<complexNum> subConstellation;
				vector<double> subProb = qam.constellationProbability(valuek, sparseH[k][k], N0, subConstellation);

				bool certain = false; int certainIdx;
				for (int s = 0; s < subProb.size(); s++)
				{
					if (subProb.at(s) == 1)
					{
						certain = true;
						certainIdx = s;
						break;
					}
				}

				if (certain == true)
				{
					path[list].at(k) = subConstellation.at(certainIdx);

					for (int j = 0; j < N; j++)
						Youtput[list].at(j) -= sparseH[j][k] * path[list].at(k);

					bufferPath.push_back(path[list]);
					bufferProb.push_back(probablity.at(list));
					bufferY.push_back(Youtput[list]);
				}
				else
				{

					for (int s = 0; s < subProb.size(); s++)
					{
						vector<complexNum> subEstimatedS(N), subYoutput(N);
						copy(path[list].begin(), path[list].end(), subEstimatedS.begin());
						copy(Youtput[list].begin(), Youtput[list].end(), subYoutput.begin());

						subEstimatedS.at(k) = subConstellation.at(s);
						for (int j = 0; j < N; j++)
							subYoutput.at(j) -= sparseH[j][k] * subConstellation.at(s);

						bufferPath.push_back(subEstimatedS);
						bufferProb.push_back(probablity.at(list)*subProb.at(s));
						bufferY.push_back(subYoutput);

					}

				}
			}

			if (bufferPath.size() > Ls)
			{
				int deleteSize = bufferPath.size() - Ls;
				for (int j = 0; j < deleteSize; j++)
				{
					int minIdx = min_element(bufferProb.begin(), bufferProb.end()) - bufferProb.begin();
					bufferPath.erase(bufferPath.begin() + minIdx, bufferPath.begin() + (minIdx + 1));
					bufferY.erase(bufferY.begin() + minIdx, bufferY.begin() + (minIdx + 1));
					bufferProb.erase(bufferProb.begin() + minIdx, bufferProb.begin() + (minIdx + 1));
				}
			}

			path.resize(bufferPath.size());
			Youtput.resize(bufferY.size());
			probablity.resize(bufferProb.size());
			for (int r = 0; r < bufferPath.size(); r++)
			{
				path[r].resize(N);
				Youtput[r].resize(N);
				copy(bufferPath[r].begin(), bufferPath[r].end(), path[r].begin());
				copy(bufferY[r].begin(), bufferY[r].end(), Youtput[r].begin());
			}
			copy(bufferProb.begin(), bufferProb.end(), probablity.begin());

			complexOperation.DeleteMatrix(Hk);

			for (int j = 0; j < N; j++)
				sparseH[j][k] = complexNum(0, 0);

			state.at(k) = 0;

			int checkSum = 0;
			for (int r = 0; r < N; r++)
				checkSum += state.at(r);
			check = checkSum;

		}
		complexOperation.DeleteMatrix(sparseH);

		int finalIndex = max_element(probablity.begin(), probablity.end()) - probablity.begin();

		vector<complexNum> hatS(N);
		copy(path[finalIndex].begin(), path[finalIndex].end(), hatS.begin());

		return hatS;
	}

private:
	complexNum PartialMMSEequalization(complexNum **Hk, vector<complexNum> Yk, double N0 ,double &channelPower)
	{
		COMPLEX_OPERATION complexOperation;

		complexNum **Her_Hk = complexOperation.HermitianOperator(Hk);
		complexNum **Mul_Hk_Her_Hk = complexOperation.MatrixMultiplication(Hk, Her_Hk);
		for (int r = 0; r <= 2 * D; r++)
			Mul_Hk_Her_Hk[r][r] += complexNum(N0, 0);
		complexNum **inv_Mul_Hk_Her_Hk = complexOperation.ComplexInverseMatrix(Mul_Hk_Her_Hk);

		vector<complexNum> Her_gk(2 * D + 1, 0 );
		double sumP = 0;
		for (int i = 0; i < Her_gk.size(); i++)
		{
			for (int j = 0; j < 2 * D + 1; j++)
				Her_gk.at(i) += Her_Hk[2 * D][j] * inv_Mul_Hk_Her_Hk[j][i];

			sumP += abs(Her_gk.at(i)*conj(Her_gk.at(i)));
		}

		channelPower = sumP;

		complexOperation.DeleteMatrix(Her_Hk);
		complexOperation.DeleteMatrix(Mul_Hk_Her_Hk);
		complexOperation.DeleteMatrix(inv_Mul_Hk_Her_Hk);

		complexNum sum = 0;
		for (int i = 0; i < Yk.size(); i++)
			sum += Her_gk.at(i)*Yk.at(i);

		return sum;
	}

	vector<unsigned int>segmentSwapIdx;

	
};

#endif // !CHANNEL_ESTIMATION

