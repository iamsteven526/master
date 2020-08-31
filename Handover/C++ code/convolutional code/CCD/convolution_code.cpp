#include"setting.h"


vector<bool> CCD::Encoder(vector<bool> message)
{
	vector<bool> codedBlock(BLOCK_LEN, 0);
	vector<bool> registerState(m + 1,0);

	for (int i = 0; i < message.size(); i++)
	{
		registerState.erase(registerState.end() - 1);
		registerState.insert(registerState.begin(), message.at(i));

		bool v1 = 0, v2 = 0;
		for (int j = 0; j < registerState.size(); j++)
		{
			v1 ^= registerState.at(j)&generatorPoly1.at(j);
			v2 ^= registerState.at(j)&generatorPoly2.at(j);
		}
		codedBlock.at(n * i) = v1;
		codedBlock.at(n*i + 1) = v2;
	}

	return codedBlock;
}

vector<double> Modulation(vector<bool> codeword) 
{
	vector<double> output(codeword.size());
	for (int i = 0; i < output.size(); i++)
		(codeword.at(i) == 0) ? output.at(i) = sqrt(Eb): output.at(i) = -sqrt(Eb);

	return output;
}

vector<double> AWGNoise(vector<double> signal ,float snr) 
{
	double R = double(k) / n;
	double linearSNR = pow(10, 0.1*snr);
	double sigma = sqrt(Eb / (2 * R *linearSNR));

	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine generator(seed);
	normal_distribution<double> normalDistribution(0, sigma);

	vector<double> receivedSignal(signal.size());
	for (int i = 0; i < receivedSignal.size(); i++)
		receivedSignal.at(i) = signal.at(i) + normalDistribution(generator);

	return receivedSignal;
}

vector<vector<double>> Quantization(vector<double> receivedSignal, vector<int> Q)
{
	vector<vector<double>> AllOutput;
	for (int q = 0; q < Q.size(); q++) 
	{
		int size = Q.at(q);
		float interval = (2*sqrt(Eb)) / size;

		vector<float> point(size);
		for (int i = 0; i < size/2; i++) 
		{
			point.at(i) = -sqrt(Eb) + i * interval;
			point.at(size - 1 - i) = -point.at(i);
		}

		vector<double> output(receivedSignal.size());
		for (int i = 0; i < receivedSignal.size(); i++)
		{
			vector<double> distance(size);
			for (int d = 0; d < size; d++)
				distance.at(d) = abs(receivedSignal.at(i) - point.at(d));

			int minIndex = min_element(distance.begin(), distance.end()) - distance.begin();
			output.at(i) = point.at(minIndex);
		}

		AllOutput.push_back(output);
	}

	return AllOutput;
}

vector<double>  CCD::ViterbiAlgorithm(vector<vector<double>> output, vector<bool> message)
{
	int stateNum = pow(2, m + 1);
	int segmentNumber = message.size();
	vector<double> bitErr;

	int Qsize = output.size();

	double par1 = 0;
	bool par2 = 0;
	for (int q = 0; q < Qsize; q++) 
	{
		int err = 0;

		double **distanceMatric = CreateMatrix(stateNum, tau, par1);		// stateNum*tau
		bool **bitPath = CreateMatrix(stateNum, tau, par2);					// stateNum*tau
		bool **tmp_bitPath = CreateMatrix(stateNum, tau, par2);

		vector<bool> decodedbit(segmentNumber, 0);

		for (int s = 0; s < segmentNumber; s++)
		{

			if (s < m + 1)
			{
				int branch = stateNum / pow(2, m + 1 - s);

				for (int j = 0; j < branch; j++)
				{
					vector<double> zeroCodeword(n);
					vector<double> oneCodeword(n);

					for (int a = 0; a < n; a++)
					{
						zeroCodeword.at(a) = zeroCodewordSet[j].at(a);
						oneCodeword.at(a) = oneCodewordSet[j].at(a);
					}

					unsigned int zeroNextState = zeroNextStateSet.at(j);
					unsigned int oneNextState = oneNextStateSet.at(j);

					if (s == 0)
					{
						distanceMatric[zeroNextState][s] = zeroCodeword.at(0)*output[q].at(2 * s) + zeroCodeword.at(1)*output[q].at(2 * s + 1);
						distanceMatric[oneNextState][s] = oneCodeword.at(0)*output[q].at(2 * s) + oneCodeword.at(1)*output[q].at(2 * s + 1);
					}
					else
					{
						for (int a = 0; a < s; a++)
						{
							bitPath[zeroNextState][a] = tmp_bitPath[j][a];
							bitPath[oneNextState][a] = tmp_bitPath[j][a];
						}
						distanceMatric[zeroNextState][s] = distanceMatric[j][s - 1] + zeroCodeword.at(0)*output[q].at(2 * s) + zeroCodeword.at(1)*output[q].at(2 * s + 1);	// 0 branch Q = 4
						distanceMatric[oneNextState][s] = distanceMatric[j][s - 1] + oneCodeword.at(0)*output[q].at(2 * s) + oneCodeword.at(1)*output[q].at(2 * s + 1);	// 1 branch Q = 4

					}

					bitPath[zeroNextState][s] = 0;
					bitPath[oneNextState][s] = 1;

				}

			}
			else
			{

				double **AugmentedDistanceMatric = CreateMatrix(2 * stateNum, tau, par1);
				bool **AugmentedBitPath = CreateMatrix(2 * stateNum, tau, par2);
				bool **AugmentedTmp_bitPath = CreateMatrix(2 * stateNum, tau, par2);
				for (int r = 0; r < stateNum; r++)
				{
					for (int c = 0; c < tau; c++)
					{
						AugmentedDistanceMatric[r][c] = distanceMatric[r][c];
						AugmentedBitPath[r][c] = bitPath[r][c];
						AugmentedTmp_bitPath[r][c] = tmp_bitPath[r][c];
					}
				}

				int loc = s;
				if (s >= tau - 1) loc = tau - 1;

				for (int j = 0; j < stateNum; j++)
				{
					vector<double> zeroCodeword(n);
					vector<double> oneCodeword(n);

					for (int a = 0; a < n; a++)
					{
						zeroCodeword.at(a) = zeroCodewordSet[j].at(a);
						oneCodeword.at(a) = oneCodewordSet[j].at(a);
					}

					unsigned int zeroNextState = zeroNextStateSet.at(j);
					unsigned int oneNextState = oneNextStateSet.at(j);

					if (j >= pow(2, m))
					{
						zeroNextState = zeroNextState + stateNum;
						oneNextState = oneNextState + stateNum;
					}

					for (int a = 0; a < loc; a++)
					{
						AugmentedBitPath[zeroNextState][a] = AugmentedTmp_bitPath[j][a];
						AugmentedBitPath[oneNextState][a] = AugmentedTmp_bitPath[j][a];
					}

					AugmentedDistanceMatric[zeroNextState][loc] = AugmentedDistanceMatric[j][loc - 1] + zeroCodeword.at(0)*output[q].at(2 * s) + zeroCodeword.at(1)*output[q].at(2 * s + 1);	// 0 branch Q = 4
					AugmentedDistanceMatric[oneNextState][loc] = AugmentedDistanceMatric[j][loc - 1] + oneCodeword.at(0)*output[q].at(2 * s) + oneCodeword.at(1)*output[q].at(2 * s + 1);	// 1 branch Q = 4

					AugmentedBitPath[zeroNextState][loc] = 0;
					AugmentedBitPath[oneNextState][loc] = 1;
				}

				for (int i = 0; i < stateNum; i++)
				{
					int choose = 0;
					if (AugmentedDistanceMatric[i][loc] >= AugmentedDistanceMatric[i + stateNum][loc])
						choose = i;
					else
						choose = i + stateNum;

					distanceMatric[i][loc] = AugmentedDistanceMatric[choose][loc];

					for (int b = 0; b < tau; b++)
						bitPath[i][b] = AugmentedBitPath[choose][b];

				}

				DeleteMatrix(AugmentedDistanceMatric);
				DeleteMatrix(AugmentedBitPath);
				DeleteMatrix(AugmentedTmp_bitPath);

				if (s >= tau - 1)
				{
					int maxIdx = 0;

					for (int i = 0; i < stateNum; i++)
					{
						if (distanceMatric[i][tau - 1] > distanceMatric[maxIdx][tau - 1])
							maxIdx = i;
					}

					if (s == segmentNumber - 1)
					{
						for (int b = 0; b < tau - 1; b++)
							decodedbit.at(s - (tau - 1) + b) = bitPath[maxIdx][b];
					}
					else
						decodedbit.at(s - (tau - 1)) = bitPath[maxIdx][0];


					for (int i = 0; i < stateNum; i++)
					{
						for (int b = 0; b < tau - 1; b++)
						{
							distanceMatric[i][b] = distanceMatric[i][b + 1];
							bitPath[i][b] = bitPath[i][b + 1];
						}
						distanceMatric[i][tau - 1] = 0;
						bitPath[i][tau - 1] = 0;

					}
				}


			}

			for (int r = 0; r < stateNum; r++)
			{
				for (int c = 0; c < tau; c++)
					tmp_bitPath[r][c] = bitPath[r][c];
			}

		}

		DeleteMatrix(distanceMatric);
		DeleteMatrix(bitPath);
		DeleteMatrix(tmp_bitPath);

		for (int b = 0; b < decodedbit.size(); b++)
		{
			if (decodedbit.at(b) != message.at(b))
				err++;
		}
		
		bitErr.push_back( err);
	}
	return bitErr;
}

