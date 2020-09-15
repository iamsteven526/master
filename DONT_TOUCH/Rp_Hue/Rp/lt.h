#pragma once

class LT
{

private:

	void LLR_I2O(int current_len);
	void LLR_O2I(int current_len);

public:

	LT(int input_len, int max_output_len, int incremental_len);

	void Encoder(int *data, int *tx);
	void Decoder(double variance, double *rx, double *app_llr, int current_len, int it_max);

	int in_len, max_out_len, ir_len, **i_edge, **o_edge, *i_deg, *o_deg, **i_order, **o_order;
	double *ch_llr, *ch_value, **llr_o2i, **llr_i2o, temp_llr[300];

};