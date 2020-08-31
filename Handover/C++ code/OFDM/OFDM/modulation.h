#ifndef MODULATION
#define MODULATION
#define _USE_MATH_DEFINES
#include"setting.h"


class PSK 
{
public:
	PSK() 
	{
		amplitude = sqrt(Es);
		constellationNum = pow(2, MOD_BIT);

		PSKTable.resize(constellationNum);
		for (int i = 1; i <= constellationNum; i++)
			PSKTable.at(i - 1) = complexNum(amplitude*cos(M_PI*(2*i-1)/ constellationNum), -amplitude*sin(M_PI*(2 * i - 1) / constellationNum));

		bool v = 0;
		garyCodingTable = CreateMatrix(constellationNum, MOD_BIT, v);
		int garyMappingBit[4] = { 0,1,1,0 };

		bitZeroSet.resize(constellationNum / 2, vector<unsigned int>(MOD_BIT));
		bitOneSet.resize(constellationNum / 2, vector<unsigned int>(MOD_BIT));

		
		for (int c = 0; c < MOD_BIT; c++) 
		{
			int size = pow(2, MOD_BIT - 1 - c); 
			int group = pow(2, c + 1);
			for (int g = 0; g < group; g++) 
			{
				for (int s = 0; s < size; s++)
					garyCodingTable[g*size + s][c] = garyMappingBit[g % 4];
			}

			int zeroIdx = 0, oneIdx = 0;
			for (int j = 0; j < constellationNum; j++)
			{
				if (garyCodingTable[j][c] == 0)
				{
					bitZeroSet[zeroIdx][c] = j;
					zeroIdx++;
				}
				else if (garyCodingTable[j][c] == 1)
				{
					bitOneSet[oneIdx][c] = j;
					oneIdx++;
				}
			}
		}
	}

	vector<int> bit2SymbolIdx(vector<int> codedbit) 
	{
		vector<int> symbolIndex(N);

		vector<int> bit2dec(constellationNum);
		for (int c = 0; c < constellationNum; c++)
		{
			int sum = 0;
			for (int b = 0; b < MOD_BIT; b++)
				sum += pow(2, b)*garyCodingTable[c][b];
			bit2dec.at(c) = sum;

		}
		
		for (int i = 0; i < N; i++) 
		{
			int codeBitDec = 0;
			for (int b = 0; b < MOD_BIT; b++)
				codeBitDec += pow(2, b)*codedbit.at(i*MOD_BIT + b);
			
			int pos = find(bit2dec.begin(), bit2dec.end(), codeBitDec)- bit2dec.begin();
			symbolIndex.at(i) = pos;
		}
		return symbolIndex;
	}

	vector<complexNum> ModulationMapping(vector<int> randomIndex)
	{
		vector<complexNum>output(N);

		for (int k = 0; k < N; k++)
			output.at(k) = PSKTable.at(randomIndex.at(k));

		return output;
	}

	int HardDemodMapping(vector<complexNum>input, vector<int> inputIndex) 
	{
		vector<int> demodMappingIndex(input.size());

		int bitDiffSum = 0;

		for (int i = 0; i < demodMappingIndex.size(); i++) 
		{
			vector<double> distance(constellationNum);
			for (int j = 0; j < constellationNum; j++)
				distance.at(j) = norm(input.at(i) - PSKTable.at(j));

			int minIdx = min_element(distance.begin(), distance.end()) - distance.begin();
			demodMappingIndex.at(i) = minIdx;


			for (int c = 0; c < MOD_BIT; c++)
			{
				if (garyCodingTable[demodMappingIndex.at(i)][c] != garyCodingTable[inputIndex.at(i)][c]) 
					bitDiffSum++;
			}
		}
		
		return bitDiffSum;
	}

	vector<double> SoftDemodMapping(vector<complexNum>input, double N0, vector<complexNum> channelFreqResponse)
	{

		vector<double> demodMappingLLR(input.size()*(MOD_BIT));
		double IQvariance = (N0 / 2);

		/*double IQvariance;
		if (N0 <= 0.0016)
			IQvariance = (0.0016 / 2);
		else
			IQvariance = (N0 / 2);*/

		for (int i = 0; i < input.size(); i++) 
		{
			vector<double> bitLLR(MOD_BIT);

			double subChannelNoiseVariance = IQvariance / norm(channelFreqResponse.at(i));

			for (int b = 0; b < MOD_BIT; b++) 
			{
				double A_dmin = 0;	//zero
				for (int c = 0; c < bitZeroSet.size(); c++) 
				{
					double distance = norm(input.at(i) - PSKTable.at(bitZeroSet.at(c).at(b)));
					if (c == 0) A_dmin = distance;
					else if (distance < A_dmin) A_dmin = distance;

				}

				double B_dmin = 0;	//zero
				for (int c = 0; c < bitOneSet.size(); c++)
				{
					double distance = norm(input.at(i) - PSKTable.at(bitOneSet.at(c).at(b)));
					if (c == 0) B_dmin = distance;
					else if (distance < B_dmin) B_dmin = distance;
				}

				bitLLR.at(b) = (B_dmin - A_dmin) / (2 * subChannelNoiseVariance);
				demodMappingLLR.at(MOD_BIT*i + b) = bitLLR.at(b);
			}
		}

		return demodMappingLLR;
	}

	vector<complexNum> PSKTable;

private:

	float amplitude;
	int constellationNum;

	bool **garyCodingTable;
	vector<vector<unsigned int>> bitZeroSet, bitOneSet;
};

class QAM
{
public:

	QAM(uint8_t bit_num) :bitNum(bit_num)
	{
		consellationNum = int(pow(2, bitNum));
		Eb = Es / bit_num;
		Eav = (3 * bitNum*Eb) / (2 * (consellationNum - 1));

		signalGain = sqrt(Eav);

		PAMTable.resize(bitNum);
		for (int i = 0; i < bitNum; i++)
			PAMTable.at(i) = ((-3)*signalGain + i * 2 * signalGain);


		threshold.resize(PAMTable.size() - 1);
		for (int i = 0; i < threshold.size(); i++)
			threshold.at(i) = (PAMTable.at(i) + PAMTable.at(i + 1)) / 2;

		PAM_bit_Zero_Aset.resize(bitNum / 2, vector<float>(bitNum / 2));
		PAM_bit_Zero_Aset.at(0).at(0) = PAMTable[0]; PAM_bit_Zero_Aset.at(0).at(1) = PAMTable[1];
		PAM_bit_Zero_Aset.at(1).at(0) = PAMTable[0]; PAM_bit_Zero_Aset.at(1).at(1) = PAMTable[3];

		PAM_bit_One_Bset.resize(bitNum / 2, vector<float>(bitNum / 2));
		PAM_bit_One_Bset.at(0).at(0) = PAMTable[2]; PAM_bit_One_Bset.at(0).at(1) = PAMTable[3];
		PAM_bit_One_Bset.at(1).at(0) = PAMTable[1]; PAM_bit_One_Bset.at(1).at(1) = PAMTable[2];

		CreateABset();
		PAM_bit_Zero_Aset.clear();	vector<vector<float>>().swap(PAM_bit_Zero_Aset);
		PAM_bit_One_Bset.clear();	vector<vector<float>>().swap(PAM_bit_One_Bset);
	}

	vector<unsigned int> bit2SymbolIdx(vector<int> codedbit)
	{
		vector<unsigned int> symbolIndex(N);

		vector<int> bit2PAMmapping((MOD_BIT / 2) * N);
		for (int i = 0; i < (MOD_BIT / 2) * N; i++)
		{

			if (codedbit.at(i*(MOD_BIT / 2)) == 0 && codedbit.at(i*(MOD_BIT / 2) + 1) == 0)
				bit2PAMmapping.at(i) = 0;
			else if (codedbit.at(i*(MOD_BIT / 2)) == 0 && codedbit.at(i*(MOD_BIT / 2) + 1) == 1)
				bit2PAMmapping.at(i) = 1;
			else if (codedbit.at(i*(MOD_BIT / 2)) == 1 && codedbit.at(i*(MOD_BIT / 2) + 1) == 1)
				bit2PAMmapping.at(i) = 2;
			else if (codedbit.at(i*(MOD_BIT / 2)) == 1 && codedbit.at(i*(MOD_BIT / 2) + 1) == 0)
				bit2PAMmapping.at(i) = 3;
		}

		for (int k = 0; k < N; k++)
		{
			int Iidx = bit2PAMmapping.at(k*(MOD_BIT / 2));
			int Qidx = bit2PAMmapping.at(k*(MOD_BIT / 2) + 1);

			symbolIndex.at(k) = Iidx * (MOD_BIT)+Qidx;
		}


		return symbolIndex;
	}

	vector<unsigned int> deMappingIdx(vector<complexNum> input)
	{
		vector<unsigned int> demodMappingIndex(input.size());

		for (int i = 0; i < input.size(); i++)
		{
			int real_idx = -1, imag_idx = -1, realBitDiff = 0, imagBitDiff = 0;
			for (int j = 0; j < threshold.size(); j++)
			{
				if (input.at(i).real() > threshold.at(j))
					real_idx = j + 1;
				else
				{
					real_idx = j;
					break;
				}
			}

			for (int j = 0; j < threshold.size(); j++)
			{
				if (input.at(i).imag() > threshold.at(j))
					imag_idx = j + 1;
				else
				{
					imag_idx = j;
					break;
				}
			}

			demodMappingIndex.at(i) = real_idx * bitNum + imag_idx;
		}

		return demodMappingIndex;
	}

	vector<complexNum> ModulationMapping(vector<unsigned int> randomIndex)
	{
		vector<complexNum>output(N);

		for (int k = 0; k < N; k++)
		{
			int Iidx = randomIndex.at(k) / bitNum;
			int Qidx = randomIndex.at(k) % bitNum;

			output.at(k) = complexNum(PAMTable.at(Iidx), PAMTable.at(Qidx));
		}

		return output;
	}

	int HardDemodMapping(vector<complexNum>intputSequence, vector<unsigned int> randomIndex)
	{
		vector<unsigned int> demodMappingIndex = deMappingIdx(intputSequence);

		int bitDiff = 0;

		for (int i = 0; i < N; i++)
		{
			int  realBitDiff = 0, imagBitDiff = 0;

			int real_part = randomIndex.at(i) / bitNum;
			int imag_part = randomIndex.at(i) % bitNum;

			int real_idx = demodMappingIndex.at(i) / bitNum;
			int imag_idx = demodMappingIndex.at(i) % bitNum;

			int real_diff = abs(real_idx - real_part);
			if (real_diff == bitNum - 1) realBitDiff = 1;
			else realBitDiff = real_diff;

			int imag_diff = abs(imag_idx - imag_part);
			if (imag_diff == bitNum - 1) imagBitDiff = 1;
			else imagBitDiff = imag_diff;

			bitDiff += (realBitDiff + imagBitDiff);

		}

		return bitDiff;
	}

	complexNum oneSampleHardDecision(complexNum input) 
	{
		vector<vector<double>> distance(2, vector<double>(PAMTable.size()));

		for (int j = 0; j < PAMTable.size(); j++)
		{
			distance[0].at(j) = norm(input.real() - PAMTable.at(j));
			distance[1].at(j) = norm(input.imag() - PAMTable.at(j));
		}

		int realIdx = min_element(distance[0].begin(), distance[0].end()) - distance[0].begin();
		int imagIdx = min_element(distance[1].begin(), distance[1].end()) - distance[1].begin();

		return complexNum(PAMTable.at(realIdx), PAMTable.at(imagIdx));
	}

	vector<double> SoftDemodMapping(vector<complexNum>intputSequence, double N0, vector<double> gainPower)
	{
		vector<double> demodMappingLLR(intputSequence.size()*(MOD_BIT));

		//double IQvariance = (N0 / 2);

		double IQvariance;
		if (N0 <= 0.0016)
			IQvariance = (0.0016 / 2);
		else
			IQvariance = (N0 / 2);

		for (int i = 0; i < intputSequence.size(); i++)
		{
			vector<double> bitLLR(MOD_BIT);

			double subChannelNoiseVariance = IQvariance * gainPower.at(i);
		
			for (int b = 0; b < bitLLR.size(); b++)
			{
				double A_dmin = 0;	//zero
				for (int c = 0; c < zero_Aset.at(0).size(); c++)
				{
					double distance = norm(intputSequence.at(i) - zero_Aset.at(b).at(c));
					if (c == 0) A_dmin = distance;
					else if (c > 0 && distance < A_dmin) A_dmin = distance;
				}

				double B_dmin = 0;	//one
				for (int c = 0; c < one_Bset.at(0).size(); c++)
				{
					double distance = norm(intputSequence.at(i) - one_Bset.at(b).at(c));
					if (c == 0) B_dmin = distance;
					else if (c > 0 && distance < B_dmin) B_dmin = distance;
				}

				bitLLR.at(b) = (B_dmin - A_dmin) / (2 * subChannelNoiseVariance);
				demodMappingLLR.at(MOD_BIT*i + b) = bitLLR.at(b);
			}
		}

		double ave = 0;
		for (auto &k : demodMappingLLR)
			ave += abs(k);

		ave = ave / demodMappingLLR.size();

		double aveGain = 0;
		for(auto &k: gainPower)
			aveGain+= k;


		//cout <<endl<< endl << "average LLR: " << ave << ", average channel Decay: " << aveGain << endl;

		return demodMappingLLR;


	}

	vector<double> constellationProbability(complexNum value, complexNum channelGain, double N0, vector<complexNum> &subConstellation)
	{
		vector<double> bitLLR(MOD_BIT);
		vector<double> bitZeroProbability(MOD_BIT), bitOneProbability(MOD_BIT);
		vector<double> InphasePAMprob_b1b2(MOD_BIT), QphasePAMprob_b3b4(MOD_BIT);
		vector<double> constellationProb(consellationNum);
		double IQvariance = (N0 / 2);
		double subChannelNoiseVariance = IQvariance / norm(channelGain);

		for (int b = 0; b < bitLLR.size(); b++)
		{
			double A_dmin = 0;	//zero
			for (int c = 0; c < zero_Aset.at(0).size(); c++)
			{
				double distance = norm(value - zero_Aset.at(b).at(c));
				if (c == 0) A_dmin = distance;
				else if (c > 0 && distance < A_dmin) A_dmin = distance;
			}

			double B_dmin = 0;	//one
			for (int c = 0; c < one_Bset.at(0).size(); c++)
			{
				double distance = norm(value - one_Bset.at(b).at(c));
				if (c == 0) B_dmin = distance;
				else if (c > 0 && distance < B_dmin) B_dmin = distance;
			}

			bitLLR.at(b) = (B_dmin - A_dmin) / (2 * subChannelNoiseVariance);
			bitOneProbability.at(b) = 1 / (exp(bitLLR.at(b)) + 1);
			bitZeroProbability.at(b) = 1 - bitOneProbability.at(b);
		}

		InphasePAMprob_b1b2.at(0) = bitZeroProbability.at(0)*bitZeroProbability.at(1);
		InphasePAMprob_b1b2.at(1) = bitZeroProbability.at(0)*bitOneProbability.at(1);
		InphasePAMprob_b1b2.at(2) = bitOneProbability.at(0)*bitOneProbability.at(1);
		InphasePAMprob_b1b2.at(3) = bitOneProbability.at(0)*bitZeroProbability.at(1);

		QphasePAMprob_b3b4.at(0) = bitZeroProbability.at(2)*bitZeroProbability.at(3);
		QphasePAMprob_b3b4.at(1) = bitZeroProbability.at(2)*bitOneProbability.at(3);
		QphasePAMprob_b3b4.at(2) = bitOneProbability.at(2)*bitOneProbability.at(3);
		QphasePAMprob_b3b4.at(3) = bitOneProbability.at(2)*bitZeroProbability.at(3);


		for (int j = 0; j < consellationNum; j++)
			constellationProb.at(j) = InphasePAMprob_b1b2.at(j / MOD_BIT)*QphasePAMprob_b3b4.at(j%MOD_BIT);

		vector<double> subProb(4);
		subConstellation.resize(4);
		cout << "value: " << value << endl;
		for (int l = 0; l < 4; l++) 
		{
			int idx = max_element(constellationProb.begin(), constellationProb.end()) - constellationProb.begin();
			subProb.at(l) = constellationProb.at(idx);
			subConstellation.at(l) = complexNum(PAMTable.at(idx / MOD_BIT), PAMTable.at(idx%MOD_BIT));

			constellationProb.at(idx) = -1;
			cout << "point: " << subConstellation.at(l) << ", prob: " << subProb.at(l) << endl;
		}

		system("pause");
		
		return subProb;
	}

	vector<complexNum> SoftMapper(vector<double> LLR, vector<double> &probability)
	{
		probability.resize(N);
		vector<complexNum> softMappedX(N);

		for (int k = 0; k < N; k++)
		{
			vector<double> bitZeroProbability(MOD_BIT), bitOneProbability(MOD_BIT);
			vector<double> InphasePAMprob_b1b2(MOD_BIT), QphasePAMprob_b3b4(MOD_BIT);

			for (int b = 0; b < MOD_BIT; b++)
			{
				bitOneProbability.at(b) = 1 / (exp(LLR.at(k*MOD_BIT + b)) + 1);
				bitZeroProbability.at(b) = 1 - bitOneProbability.at(b);
			}

			InphasePAMprob_b1b2.at(0) = bitZeroProbability.at(0)*bitZeroProbability.at(1);
			InphasePAMprob_b1b2.at(1) = bitZeroProbability.at(0)*bitOneProbability.at(1);
			InphasePAMprob_b1b2.at(2) = bitOneProbability.at(0)*bitOneProbability.at(1);
			InphasePAMprob_b1b2.at(3) = bitOneProbability.at(0)*bitZeroProbability.at(1);

			QphasePAMprob_b3b4.at(0) = bitZeroProbability.at(2)*bitZeroProbability.at(3);
			QphasePAMprob_b3b4.at(1) = bitZeroProbability.at(2)*bitOneProbability.at(3);
			QphasePAMprob_b3b4.at(2) = bitOneProbability.at(2)*bitOneProbability.at(3);
			QphasePAMprob_b3b4.at(3) = bitOneProbability.at(2)*bitZeroProbability.at(3);

			vector<double> constellationProb(consellationNum);
			complexNum mappedX = 0;
			for (int j = 0; j < consellationNum; j++)
			{
				constellationProb.at(j) = InphasePAMprob_b1b2.at(j / MOD_BIT)*QphasePAMprob_b3b4.at(j%MOD_BIT);
				mappedX += constellationProb.at(j)*complexNum(PAMTable.at(j / MOD_BIT), PAMTable.at(j%MOD_BIT));
			}
			int idx = max_element(constellationProb.begin(), constellationProb.end()) - constellationProb.begin();
			probability.at(k) = constellationProb.at(idx);
			softMappedX.at(k) = mappedX;
		}


		return softMappedX;
	}

	vector<float> PAMTable;
	float signalGain;

private:
	uint8_t bitNum;
	float Eb, Eav;
	int consellationNum;

	vector<float> threshold;
	vector<vector<float>> PAM_bit_Zero_Aset;	// row0: b0 = 0  row2:b1 = 0;
	vector<vector<float>> PAM_bit_One_Bset;		// row0: b0 = 1  row2:b1 = 1;

	vector<vector<complexNum>> zero_Aset;
	vector<vector<complexNum>> one_Bset;

	void CreateABset()
	{
		zero_Aset.resize(0, vector<complexNum>(0));				// row: 4 bit col :  8 corresponding constellation points
		one_Bset.assign(0, vector<complexNum>(0));

		for (int b = 0; b < MOD_BIT; b++)
		{
			vector<complexNum> onebit_zero_Aset(0);
			vector<complexNum> onebit_one_Bset(0);

			if (b < MOD_BIT / 2)	// b0,b1: Inphase
			{
				for (int n = 0; n < MOD_BIT / 2; n++)
				{
					for (int k = 0; k < MOD_BIT; k++)
						onebit_zero_Aset.push_back(complexNum(PAM_bit_Zero_Aset.at(b).at(n), PAMTable.at(k)));
				}

				for (int n = 0; n < MOD_BIT / 2; n++)
				{
					for (int k = 0; k < MOD_BIT; k++)
						onebit_one_Bset.push_back(complexNum(PAM_bit_One_Bset.at(b).at(n), PAMTable.at(k)));
				}
			}
			else  // b2,b3:  Quadrature
			{
				for (int n = 0; n < MOD_BIT / 2; n++)
				{
					for (int k = 0; k < MOD_BIT; k++)
						onebit_zero_Aset.push_back(complexNum(PAMTable.at(k), PAM_bit_Zero_Aset.at(b - (MOD_BIT / 2)).at(n)));
				}

				for (int n = 0; n < MOD_BIT / 2; n++)
				{
					for (int k = 0; k < MOD_BIT; k++)
						onebit_one_Bset.push_back(complexNum(PAMTable.at(k), PAM_bit_One_Bset.at(b - (MOD_BIT / 2)).at(n)));
				}
			}
			zero_Aset.push_back(onebit_zero_Aset);
			one_Bset.push_back(onebit_one_Bset);
		}
	}

};



#endif // !MODULATION

