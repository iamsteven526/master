#include"setting.h"
#include"math_operation.h"
#include"modulation.h"
#include"channel_coding.h"
#include"channel_model.h"
#include"channel_estimation.h"

ESTIMATION channelEstimation;
LDPC ldpc;
CHANNEL channel;

void FileSetting();
vector<complexNum> Transmitter(vector<unsigned int> &randomIndex);
vector<complexNum> Channel(vector<complexNum> TtransmitterIFFTSignal, vector<vector<complexNum>> cyclicPrefix, vector<vector<complexNum>> &nextSymbolOverhead, vector<complexNum> &channelImpulseResponse, vector<vector<complexNum>> &timeVaryingChannelImpulseResponse, vector<double> &oneCFRlPower, vector<double> &fastFadingPower, double N0, int OneFrameSymbol);
double Receiver(vector<complexNum> channelImpulseResponse, vector<complexNum> frequencyDomainReceivedSignal, vector<vector<complexNum>>timeVaryingChannelImpulseResponse, vector<unsigned int> randomIndex, vector<double> &flagNum, double N0);

int main() 
{
	FileSetting();

	vector<double> oneCFRPower, fastFadingPower;
	vector<complexNum> channelImpulseResponse;

	int SNR_point = (SNR_dB_END - SNR_dB_START) / INCREMENT + 1, snr_point = 0;
	vector<double> bler(SNR_point), ber(SNR_point);
	vector<double> LS_bler_0th(SNR_point), LS_bler_1st(SNR_point), LS_bler_2nd(SNR_point), LS_ber_0th(SNR_point), LS_ber_1st(SNR_point), LS_ber_2nd(SNR_point);
	vector<double> LS_dia_mse(SNR_point), LS_matrix_mse_1st(SNR_point), LS_matrix_mse_2nd(SNR_point);
	float snr = SNR_dB_START;

	for (; snr <= SNR_dB_END; snr += INCREMENT) 
	{
		double linearSNR = pow(10, 0.1*snr), N0 = Es / (MOD_BIT*R*linearSNR);			
		cout << endl << endl << "SNR: " << snr << "(dB), variance N0: " << N0 << endl;

		/* initialization*/
		srand(time(NULL));
		vector<unsigned int> randomIndex;
		vector<double> softLLR, flagNum(4, 0);
		vector<complexNum> timeDomainTransmittedSignal, diagonalElement, frequencyDomainReceivedSignal;
		vector<vector<complexNum>>timeVaryingChannelImpulseResponse, cyclicPrefix(2, vector<complexNum>(G)), nextSymbolOverhead(2, vector<complexNum>(L - 1, 0));
		unsigned int symbolNum = 0, OneFrameSymbol = 0, allBlockError = 0, allBitError = 0;
		
		while (allBitError < MIN_BIT_ERR || allBlockError < MIN_BLK_ERR || symbolNum / SymbolPerFrame < RadioFrameNum)
		{
			symbolNum++;
			OneFrameSymbol = symbolNum % SymbolPerFrame;

			timeDomainTransmittedSignal = Transmitter(randomIndex);

			frequencyDomainReceivedSignal = Channel(timeDomainTransmittedSignal, cyclicPrefix, nextSymbolOverhead, channelImpulseResponse, timeVaryingChannelImpulseResponse, oneCFRPower, fastFadingPower, N0, OneFrameSymbol);

			int bitError = Receiver(channelImpulseResponse, frequencyDomainReceivedSignal, timeVaryingChannelImpulseResponse, randomIndex, flagNum, N0);
			
			if (bitError) allBlockError++;
			allBitError += bitError;
			double BLER = double(allBlockError) / symbolNum;
			double BER = double(allBitError) / (symbolNum*N*MOD_BIT*R);

			double averagePower = 0;
			for (int l = 0; l < oneCFRPower.size(); l++)
				averagePower += oneCFRPower.at(l);
			averagePower = averagePower / oneCFRPower.size();

			double fastAveragePower = 0;
			for (int l = 0; l < fastFadingPower.size(); l++)
				fastAveragePower += fastFadingPower.at(l);
			fastAveragePower = fastAveragePower / fastFadingPower.size();

			copy(cyclicPrefix[1].begin(), cyclicPrefix[1].end(), cyclicPrefix[0].begin());
			copy(nextSymbolOverhead[1].begin(), nextSymbolOverhead[1].end(), nextSymbolOverhead[0].begin());

			cout << "\rsymbol#" << symbolNum << ", Pav: " << averagePower << ", (ideal) BLER: " << BLER <<", BER: " << BER;
		}

		for (int i = 0; i < flagNum.size(); i++)
			flagNum.at(i) /= symbolNum;

		double BLER = double(allBlockError) / symbolNum;
		double BER = double(allBitError) / (symbolNum*N * MOD_BIT*R);
		bler.at(snr_point) = BLER;			ber.at(snr_point) = BER;

		FILE *subresult;
		subresult = fopen("result.txt", "a");
		fprintf(subresult, "SNR: %.1f(dB)\n", snr);
		fprintf(subresult, "Symbol num: %d, (ideal)BLER: %.9f, BER: %.9f\n", symbolNum, BLER, BER);
		fprintf(subresult, "(flag prob.) flag0: %.9f, flag1: %.9f, flag2: %.9f\n", flagNum.at(0), flagNum.at(1), flagNum.at(2));
		fclose(subresult);
		snr_point++;
	}

	//system("pause");
	return 0;
}

void FileSetting() 
{
	vector<string> channelList;
	channelList.push_back("AWGN channel");
	channelList.push_back("Slow fading- exponentially decay channel");
	channelList.push_back("Fast fading- time-varying Rayleigh fading channel");

	cout << pow(2, MOD_BIT) << "PSK-OFDM, LDPC coded" << endl;
	cout << "One-tap equalizer" << endl << endl;
	cout << "1). subcarrier number: " << N << ", subcarrier spacing: " << DELRA_FREQ / 1000 << "(kHz)" << endl;
	cout << "2). path number : " << L << "(tapes), UE speed: " << UE_SPEED << "(km/h)" << ", doppler freq: " << CARRIER_FREQ * (UE_SPEED*(1000. / 3600)) / (3. * pow(10, 8))*SYMBOL_DURATION << endl;
	cout << "3). symbol duration = " << SYMBOL_DURATION << endl;
	cout << "4). SNR_dB range: " << SNR_dB_START << "-" << SNR_dB_END << "(dB)" << endl;
	cout << "5). Channel type: " << channelList.at(CHANNEL_TYPE) << endl;
	cout << "6). Chip sample rate: " << float(Fs) / 1000000 << "(MHz)" << endl;

	FILE *result = fopen("result.txt", "w");
	fprintf(result, "OFDM system %dPSK\n", pow(2, MOD_BIT));
	fprintf(result, "One-tap equalizer  \n\n");
	fprintf(result, "1). subcarrier number: %d, subcarrier spacing: %d(kHz) \n", N, DELRA_FREQ / 1000);
	fprintf(result, "2). path number : %d(tapes), UE speed: %d(km/h)\n", L, UE_SPEED);
	fprintf(result, "3). %s\n\n", channelList.at(CHANNEL_TYPE).c_str());
	fclose(result);
}

vector<complexNum> Transmitter(vector<unsigned int> &randomIndex)
{
	vector<complexNum> frequencyDomainTransmittedSignal, timeDomainTransmittedSignal;

	if (R != 1)
	{

		vector<complexNum> modulatedSignal = qam.ModulationMapping(qam.bit2SymbolIdx(ldpc.MutlipleSubBlockEncoder()));
		vector<complexNum> frequencyDomainTransmittedSignal = channelEstimation.Permutation(modulatedSignal);
		randomIndex = qam.deMappingIdx(frequencyDomainTransmittedSignal);
	}
	else 	//uncode
	{
		randomIndex.resize(N);																			
		for (int k = 0; k < N; k++)
			randomIndex.at(k) = rand() % int(pow(2, MOD_BIT));

		frequencyDomainTransmittedSignal = qam.ModulationMapping(randomIndex);
	}
	timeDomainTransmittedSignal = fft.NormalizedInverseFastFourierTransform(frequencyDomainTransmittedSignal);

	return timeDomainTransmittedSignal;
}

vector<complexNum> Channel(vector<complexNum> TtransmitterIFFTSignal, vector<vector<complexNum>> cyclicPrefix, vector<vector<complexNum>> &nextSymbolOverhead,vector<complexNum> &channelImpulseResponse, vector<vector<complexNum>> &timeVaryingChannelImpulseResponse, vector<double> &oneCFRlPower, vector<double> &fastFadingPower,double N0, int OneFrameSymbol)
{
	vector<complexNum>receivedSignalCPRemoval, TtransmitterIFFTSignalWithCP, receivedSignal, currentSymbolwithCP, complexGaussianNoise;

	switch (CHANNEL_TYPE) {
	case 0:
	{
		receivedSignalCPRemoval.resize(N);
		for (int i = 0; i < N; i++)
			receivedSignalCPRemoval.at(i) = TtransmitterIFFTSignal.at(i);

	}break;
	case 1:
	case 2:
	{
		copy(TtransmitterIFFTSignal.begin() + (N - G), TtransmitterIFFTSignal.end(), cyclicPrefix[1].begin());

		TtransmitterIFFTSignalWithCP.resize(N);
		copy(TtransmitterIFFTSignal.begin(), TtransmitterIFFTSignal.end(), TtransmitterIFFTSignalWithCP.begin());
		TtransmitterIFFTSignalWithCP.insert(TtransmitterIFFTSignalWithCP.begin(), cyclicPrefix[1].begin(), cyclicPrefix[1].end());  // note: inserting will increase the length of vector (critical)

		if (CHANNEL_TYPE == 1)
		{
			if (OneFrameSymbol == 1)
				channelImpulseResponse = channel.ExponentiallyDecayChannel(L, ALPHA);
			receivedSignal = channel.QuasiStaticChannel(TtransmitterIFFTSignalWithCP, channelImpulseResponse, oneCFRlPower, OneFrameSymbol);
		}
		else if (CHANNEL_TYPE == 2)
		{
			if (OneFrameSymbol == 1)
			{
				channelImpulseResponse = channel.ExponentiallyDecayChannel(L, ALPHA);
				fastFadingPower.resize(0);
			}
			timeVaryingChannelImpulseResponse = channel.TimeVaryingRayleighFadingChannel(channelImpulseResponse, OneFrameSymbol);
			receivedSignal = channel.FastFadingChannel(TtransmitterIFFTSignalWithCP, timeVaryingChannelImpulseResponse, fastFadingPower);
			
			if (OneFrameSymbol == 0)
			{
				double sum = 0;
				for (int i = 0; i < fastFadingPower.size(); i++)
					sum += fastFadingPower.at(i);
				oneCFRlPower.push_back(sum / fastFadingPower.size());
			}
		}

		currentSymbolwithCP.resize(N + G);
		copy(receivedSignal.begin(), receivedSignal.begin() + (N + G), currentSymbolwithCP.begin());

		copy(receivedSignal.begin() + (N + G), receivedSignal.end(), nextSymbolOverhead[1].begin());


		for (int i = 0; i < nextSymbolOverhead[0].size(); i++)
			currentSymbolwithCP.at(i) += nextSymbolOverhead[0].at(i);

		receivedSignalCPRemoval.resize(N);
		copy(currentSymbolwithCP.begin() + G, currentSymbolwithCP.end(), receivedSignalCPRemoval.begin());

	}break;
	}

	complexGaussianNoise = channel.AWGNNoise(N, N0);
	for (int i = 0; i < receivedSignalCPRemoval.size(); i++)
		receivedSignalCPRemoval.at(i) += complexGaussianNoise.at(i);

	 vector<complexNum> outputVector = fft.NormalizedFastFourierTransform(receivedSignalCPRemoval);

	 return outputVector;
}

double Receiver(vector<complexNum> channelImpulseResponse, vector<complexNum> frequencyDomainReceivedSignal, vector<vector<complexNum>>timeVaryingChannelImpulseResponse, vector<unsigned int> randomIndex,vector<double> &flagNum, double N0)
{
	int bitErr = 0;
	vector<complexNum> CFR = fft.FastFourierTransform(channelImpulseResponse);
	vector<complexNum> estimatedOutput; 
	vector<double> processGainPower, permutedProcessGainPower;
	switch (CHANNEL_TYPE) {
	case 0: {
		estimatedOutput.resize(N);
		copy(frequencyDomainReceivedSignal.begin(), frequencyDomainReceivedSignal.end(), estimatedOutput.begin());

	}break;
	case 1:
		estimatedOutput = channelEstimation.OneTapEqualizer(CFR, frequencyDomainReceivedSignal);
		break;
	case 2:
	{
		complexNum  **H = channelEstimation.CreateHmatrix(timeVaryingChannelImpulseResponse, G);
		estimatedOutput = channelEstimation.PartialMMSEequalizer(H, frequencyDomainReceivedSignal, N0, processGainPower);
		complexOperation.DeleteMatrix(H);
	}break;
	}

	if (R != 1)
	{
		ldpc.checkTable.resize(0);	ldpc.errTable.resize(0); ldpc.errFlag.resize(0);
		vector<double> idealDecodedLLR;
		vector<complexNum> InversePermutated_X = channelEstimation.InversePermutation(estimatedOutput, processGainPower, permutedProcessGainPower);
		vector<double> softLLR = qam.SoftDemodMapping(InversePermutated_X, N0, permutedProcessGainPower);
		vector<int> decodeCodeword = ldpc.MutlipleSubBlockDecodeer(softLLR, idealDecodedLLR, bitErr, 0);

		switch (ldpc.errFlag.at(0))
		{
		case 0:
			flagNum.at(0)++;	break;
		case DEVISION:
			flagNum.at(2)++;	break;
		default:
			flagNum.at(1)++;	break;
		}
	}
	else
		bitErr = qam.HardDemodMapping(estimatedOutput, randomIndex);

	return bitErr;
}