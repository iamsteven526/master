#include <iostream>
#include "ldpc.h"
#include "lt.h"
#include "parameters.h"
using namespace std;

#pragma warning(disable:4996)

int main()
{
	LDPC ldpc(LDPC_H_COL, LDPC_H_ROW);
	LT lt(INTER_CODE_LEN, MAX_CODE_LEN, IR_LEN);
	FILE *result_txt = fopen("result.txt", "w");
	int *data = new int[DATA_LEN];
	int *interSeq = new int[INTER_CODE_LEN];
	int *tx = new int[MAX_CODE_LEN];
	double *rx = new double[MAX_CODE_LEN];
	double *interLlr = new double[INTER_CODE_LEN];
	double *decodedLlr = new double[INTER_CODE_LEN];
	double tPut[SNR_NUM], eff[SNR_NUM];
	printf("Raptor code\n\n");
	fprintf(result_txt, "Raptor code\n\n");
	if (SYSTEMATIC == 0)
	{
		printf("Non-systematic, IR length: %d\n", IR_LEN);
		fprintf(result_txt, "Non-systematic, IR length: %d\n", IR_LEN);
	}
	if (SYSTEMATIC == 1)
	{
		printf("Systematic, IR length: %d\n", IR_LEN);
		fprintf(result_txt, "Systematic, IR length: %d\n", IR_LEN);
	}
	printf("degree set: %d\n", DEGREE_SET);
	fprintf(result_txt, "degree set: %d\n", DEGREE_SET);
	if (DEC_SCHEME == 1)
	{
		printf("tandem decoding, iteration#: LDPC %d, LT %d\n\n", LDPC_IT, LT_IT);
		fprintf(result_txt, "tandem decoding, iteration#: LDPC %d, LT %d\n\n", LDPC_IT, LT_IT);
	}
	if (DEC_SCHEME == 2)
	{
		printf("joint decoding, iteration#: %d\n\n", JOINT_IT);
		fprintf(result_txt, "joint decoding, iteration#: %d\n\n", JOINT_IT);
	}
	for (int i = 0; i < SNR_NUM; i++)
	{
		double snrdB = SNR_START + (double)(SNR_STEP*i); // SNR per transmitted symbol
		double variance = 0.5*pow(10, -(snrdB / 10));
		double capacity = Capacity(variance);
		int ack, irCount, totalMsg = 0, reqSymbol = 0;
		printf("SNR[dB] = %.1f, Capacity = %f\n", snrdB, capacity);
		for (int block = 1; block <= BLOCK_NUM; block++)
		{
			RaptorEncoder(ldpc, lt, data, interSeq, tx);
			Channel(sqrt(variance), tx, rx);
			if (DEC_SCHEME == 1) TandemDecoder(lt, ldpc, variance, rx, interLlr, decodedLlr, data, ack, irCount);
			if (DEC_SCHEME == 2) JointDecoder(lt, ldpc, variance, rx, decodedLlr, data, ack, irCount);
			if (ack == 1) totalMsg += DATA_LEN;
			reqSymbol += INTER_CODE_LEN + irCount * IR_LEN;
			
			
			tPut[i] = totalMsg / (double)reqSymbol;
			eff[i] = tPut[i] / capacity;
			printf("Block#= %d, Tput = %f, Eff. = %f\r", block, tPut[i], eff[i]);
		}
		fprintf(result_txt, "SNR[dB] = %.1f, Capacity = %f, Tput = %f, Eff. = %f\n", snrdB, capacity, tPut[i], eff[i]);
		printf("\n\n");
	}
	fprintf(result_txt, "\n");
	for (int i = 0; i < SNR_NUM; i++)
	{
		fprintf(result_txt, "%f, ", tPut[i]);
	}
	fprintf(result_txt, "\n\n");
	for (int i = 0; i < SNR_NUM; i++)
	{
		fprintf(result_txt, "%f, ", eff[i]);
	}
	fprintf(result_txt, "\n");
	fclose(result_txt);
	system("pause");
	return 0;
}