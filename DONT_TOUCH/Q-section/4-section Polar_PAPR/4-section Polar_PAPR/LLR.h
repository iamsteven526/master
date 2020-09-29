

void LLR(double Recovered_Mod_real, double Recovered_Mod_imag, double sigmasquare, double LLR[p]) {

	int bit_idx;

	LLR[0] = (exp(-pow(Recovered_Mod_real + 3.0 / sqrt(10), 2) / (2 * sigmasquare)) + exp(-pow(Recovered_Mod_real + 1.0 / sqrt(10), 2) / (2 * sigmasquare))) / (exp(-pow(Recovered_Mod_real - 1.0 / sqrt(10), 2) / (2 * sigmasquare)) + exp(-pow(Recovered_Mod_real - 3.0 / sqrt(10), 2) / (2 * sigmasquare)));

	LLR[1] = (exp(-pow(Recovered_Mod_real + 3.0 / sqrt(10), 2) / (2 * sigmasquare)) + exp(-pow(Recovered_Mod_real - 3.0 / sqrt(10), 2) / (2 * sigmasquare))) / (exp(-pow(Recovered_Mod_real + 1.0 / sqrt(10), 2) / (2 * sigmasquare)) + exp(-pow(Recovered_Mod_real - 1.0 / sqrt(10), 2) / (2 * sigmasquare)));

	LLR[2] = (exp(-pow(Recovered_Mod_imag - 3.0 / sqrt(10), 2) / (2 * sigmasquare)) + exp(-pow(Recovered_Mod_imag - 1.0 / sqrt(10), 2) / (2 * sigmasquare))) / (exp(-pow(Recovered_Mod_imag + 1.0 / sqrt(10), 2) / (2 * sigmasquare)) + exp(-pow(Recovered_Mod_imag + 3.0 / sqrt(10), 2) / (2 * sigmasquare)));

	LLR[3] = (exp(-pow(Recovered_Mod_imag - 3.0 / sqrt(10), 2) / (2 * sigmasquare)) + exp(-pow(Recovered_Mod_imag + 3.0 / sqrt(10), 2) / (2 * sigmasquare))) / (exp(-pow(Recovered_Mod_imag - 1.0 / sqrt(10), 2) / (2 * sigmasquare)) + exp(-pow(Recovered_Mod_imag + 1.0 / sqrt(10), 2) / (2 * sigmasquare)));

	for (bit_idx = 0; bit_idx < p; bit_idx++) {
		LLR[bit_idx] = log(LLR[bit_idx]);
		if (LLR[bit_idx] == 0) LLR[bit_idx] = 0.000000001;
	}

}

