#include<iostream>
#include<vector>
#include<iomanip>
#include"parameter.h"
#include <algorithm>
using namespace std;


int main()
{

	double variance;
	int recursive = 0;
	bool triger = true;
	int error=0;
	int reg = 0;

	if (DEC_SCHEME == 1)
		cout << "Raptor code"<<endl<<"tandem decoding" << endl;
	else
		cout << "Raptor code" << endl<<"joint decoding" << endl;


	vector<double>rx(MAX_CODE_LEN, 0);
	vector<double>decoded_rx(INTER_CODE_LEN,0);

	LDPC ldpc(LDPC_H_COL,LDPC_H_ROW);
	LT lt(INTER_CODE_LEN, MAX_CODE_LEN);
	Raptor raptor(ldpc, lt, DATA_LEN, INTER_CODE_LEN, MAX_CODE_LEN);
	

	vector<int> random_message(DATA_LEN);
	generate(random_message.begin(), random_message.end(), [] { return rand()%2; });

	for (double snrdB = SNR_START; snrdB  < SNR_START + SNR_NUM * SNR_STEP; snrdB += SNR_STEP)
	{
		cout << endl<< "--------------------------" << endl<<"SNR[db]=" << snrdB << endl << "--------------------------" << endl;
		variance = 0.5*pow(10, -(snrdB / 10));

		for (int i = 1; i <= BLOCK_NUM; i++)
		{
			if (i % 20 == 0)
				generate(random_message.begin(), random_message.end(), [] { return rand() % 2; });

			vector<int> tx_ldpc = ldpc.encoder(random_message);
			for (recursive; recursive < IR_MAX; recursive++)
			{
				error = 0;
				vector<int> tx_lt = lt.encoder(tx_ldpc, recursive);
				channel(sqrt(variance), tx_lt, rx, recursive);
				if (DEC_SCHEME == 1)
					decoded_rx = raptor.tandem_decoder(snrdB, rx, recursive, triger, ldpc, lt);		//tandem decoder
				else
					decoded_rx = raptor.joint_decoder(recursive, snrdB, rx, triger, ldpc);	//joint decoder

				if (triger == 0)
				{
					for (int i = 0; i < DATA_LEN; i++)
						if (HARD(decoded_rx[i + INTER_CODE_LEN - DATA_LEN]) != random_message[i])
							error++;
					if (error == 0)
						break;
				}
				//if (i % 20 == 0)
				//	cout << (double)i/BLOCK_NUM<<"%" << " " << recursive << "\r";
			}
			reg += recursive;
		//	if(i%20==0)
				cout << (double)i / BLOCK_NUM*100 << "%"<<"   Tput=" << DATA_LEN * i / (INTER_CODE_LEN*i + double(reg)*IR_LEN) << "\r";
			recursive = 0;
		}
		reg = 0;

	}


	rx.clear();
	decoded_rx.clear();
	random_message.clear();
	cout << endl;

	system("pause");
	return 0;
}
