#pragma once

class LDPC
{

private:

	void Swap(unsigned short **x, unsigned short **y) 
	{
		unsigned short *temp;
		temp = *x;
		*x = *y;
		*y = temp;
	}

public:

	LDPC(int column, int row);

	void Encoder(int *data, int *tx);
	void Decoder(double variance, double *rx, double *decoded_llr, int it_max);

	unsigned short **G, *col_w, *row_w, **col_edge, **row_edge;
	int h_col, h_row, data_len, code_len;
	double *ch_llr, **app_llr, **ext_llr, temp_llr[300];

};