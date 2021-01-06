#pragma once

class LDPC
{

private:

	void swap(unsigned short **x, unsigned short **y) 
	{
		unsigned short *temp;
		temp = *x;
		*x = *y;
		*y = temp;
	}

	void DiffGraphConstruction();
	void JointI2O(double *input, double *output);
	void JointO2I();

	int out_len, *i_deg, *o_deg, **i_edge, **o_edge, **i_link, **o_link;
	double *ch_value, **llr_o2i, **llr_i2o, *input_llr, *output_llr;
	double ***p_i2o, ***p_o2i, **input_p, **output_p;

public:

	LDPC(int column, int row);

	void Encoder(int *data, int *tx);

	void SPA(double *rx, double *decoded_llr, int it_max);
	void GSPA(double **app, int it);

	void JointSPA(double *appLlr, double *refLlr, double *decodedLlr, int outer_it_max, int inter_it_max);
	void JointGSPA(double **app, int outer_it_max, int inter_it_max);

	unsigned short **G, *col_w, *row_w, **col_edge, **row_edge;
	int h_col, h_row, data_len, code_len;
	double ***p_v2c, ***p_c2v;
	double *ch_llr, **app_llr, **ext_llr, temp_value[300];

};