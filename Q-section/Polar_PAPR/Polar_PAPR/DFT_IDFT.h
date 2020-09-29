

///////////////IDFT///////////////
void IDFT(double x_real[L * N], double x_imag[L * N], double y_real[L * N], double y_imag[L * N]) {


	int		x_idx, y_idx;


	for (y_idx = 0; y_idx < L * N; y_idx++) {
		y_real[y_idx] = 0;
		y_imag[y_idx] = 0;



		for (x_idx = 0; x_idx < L * N; x_idx++) {
			y_real[y_idx] += x_real[x_idx] * cos(2.0 * PI * x_idx * y_idx / (L * N)) - x_imag[x_idx] * sin(2.0 * PI * x_idx * y_idx / (L * N));
			y_imag[y_idx] += x_real[x_idx] * sin(2.0 * PI * x_idx * y_idx / (L * N)) + x_imag[x_idx] * cos(2.0 * PI * x_idx * y_idx / (L * N));
		}
		y_real[y_idx] = (1.0 / sqrt(L * N)) * y_real[y_idx];
		y_imag[y_idx] = (1.0 / sqrt(L * N)) * y_imag[y_idx];

	}


}



///////////////DFT///////////////
void DFT(double x_real[L * N], double x_imag[L * N], double y_real[L * N], double y_imag[L * N]) {

	int		x_idx, y_idx;


	for (y_idx = 0; y_idx < L * N; y_idx++) {
		y_real[y_idx] = 0;
		y_imag[y_idx] = 0;




		for (x_idx = 0; x_idx < L * N; x_idx++) {
			y_real[y_idx] += x_real[x_idx] * cos(2.0 * PI * y_idx * x_idx / (L * N)) + x_imag[x_idx] * sin(2.0 * PI * y_idx * x_idx / (L * N));
			y_imag[y_idx] += x_real[x_idx] * -sin(2.0 * PI * y_idx * x_idx / (L * N)) + x_imag[x_idx] * cos(2.0 * PI * y_idx * x_idx / (L * N));
		}
		y_real[y_idx] = (1.0 / sqrt(L * N)) * y_real[y_idx];
		y_imag[y_idx] = (1.0 / sqrt(L * N)) * y_imag[y_idx];


	}


}




