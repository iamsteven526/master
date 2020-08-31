#include"setting.h"

int main() 
{
	CCD ccd;
	vector<int> Qsize;
	Qsize.push_back(2);	Qsize.push_back(4);	Qsize.push_back(8);

	int SNR_point = (SNR_dB_END - SNR_dB_START) / INCREMENT + 1;
	
	FILE *result = fopen("result(b).txt", "w");
	fprintf(result, "(%d, %d, %d) convolution code with Q =2/4/8 and truncation length =%d \n\n", n, k, m, tau);
	fclose(result);

	cout << "(" << n << ", " << k << ", " << m << ") convolution code with Q =2/4/8 and truncation length = " << tau << endl << endl;
	for (float snr = SNR_dB_START; snr <= SNR_dB_END; snr += INCREMENT)
	{
		cout << endl << "SNR: " << snr << "(dB)" << endl;

		srand(time(NULL));
		vector<bool> message(L + m, 0);
		vector<bool> codeword;
		vector<double> signal;
		vector<double>receivedSignal; 
		vector<vector<double>> output;
		vector<double> singleBitErr(Qsize.size());
		vector<double> sumBitErr(Qsize.size(), 0);

		unsigned int block = 0, blockErr = 0;
		while (blockErr < MIN_BLOCK_ERR || sumBitErr.at(2) < MIN_BIT_ERR)
		{
			block++;

			for (int i = 0; i < L; i++)
				message.at(i) = rand() % 2;

			codeword = ccd.Encoder(message);
			signal = Modulation(codeword);

			receivedSignal = AWGNoise(signal, snr);

			output = Quantization(receivedSignal, Qsize);
			singleBitErr = ccd.ViterbiAlgorithm(output, message);

			for (int q = 0; q < sumBitErr.size(); q++)
				sumBitErr.at(q) += singleBitErr.at(q);

			if (sumBitErr.at(2) != 0)blockErr++;
			
			cout << "\rblock#" << block << ", Q size = 2/4/8, BER = " << sumBitErr.at(0) / (INFO_LEN * block) << "/" << sumBitErr.at(1) / (INFO_LEN * block) << "/" << sumBitErr.at(2) / (INFO_LEN * block);
		}
		for (int q = 0; q < sumBitErr.size(); q++)
			sumBitErr.at(q) = sumBitErr.at(q) / (INFO_LEN * block);

	

		FILE *subresult;
		subresult = fopen("result(b).txt", "a");
		fprintf(subresult, "SNR: %.1f(dB), block#%d, BER: %.9f/%.9f/%.9f\n", snr, block, sumBitErr.at(0), sumBitErr.at(1), sumBitErr.at(2));
		fclose(subresult);
	}



	system("pause");
	return 0;

}

