
using namespace std;


void ACE_SGP(double S_real[L * N], double S_imag[L * N], double s_real[L * N], double s_imag[L * N]) {

	int			OFDM_sym_idx, mod_sym_idx, iter_idx;
	double		power[L * N] = { 0 };
	double		avg_pow = 0;
	double		A;									// clip_threshold = 4 dB; 


	double		amplitude[L * N] = { 0 };
	double		sita = 0;
	double		phase = 0;

	double		s_clip_real[L * N] = { 0 };
	double		s_clip_imag[L * N] = { 0 };

	double		c_real[L * N] = { 0 };
	double		c_imag[L * N] = { 0 };
	double		C_real[L * N] = { 0 };
	double		C_imag[L * N] = { 0 };


	double		E;
	int			n_max;
	double		c_proj[L * N] = { 0 };
	double		delta[L * N] = { 0 };
	double		Mu = 0;









	///////////////OFDM Signal Modulation///////////////
	IDFT(S_real, S_imag, s_real, s_imag);


	///////////////determining the clipping threshold A///////////////
	for (iter_idx = 0; iter_idx < 2; iter_idx++) {
		for (OFDM_sym_idx = 0; OFDM_sym_idx < L * N; OFDM_sym_idx++) {
			power[OFDM_sym_idx] = pow(s_real[OFDM_sym_idx], 2) + pow(s_imag[OFDM_sym_idx], 2);
			avg_pow += power[OFDM_sym_idx];
		}
		avg_pow /= (double)L * N;
		//cout << avg_pow << endl; system("pause");
		A = avg_pow * pow(10, 0.4);					//clipping threshold = 4dB


		E = 0;
		for (OFDM_sym_idx = 0; OFDM_sym_idx < L * N; OFDM_sym_idx++) {

			amplitude[OFDM_sym_idx] = sqrt(power[OFDM_sym_idx]);
			//cout << amplitude[OFDM_sym_idx] << endl; system("pause");
			if (amplitude[OFDM_sym_idx] > E) {
				E = amplitude[OFDM_sym_idx];
				n_max = OFDM_sym_idx;
			}
			sita = atan(s_imag[OFDM_sym_idx] / s_real[OFDM_sym_idx]);
			if (s_real[OFDM_sym_idx] < 0) phase = sita + PI;
			else phase = sita;
			//cout << amplitude[OFDM_sym_idx] << " "; system("pause");


			///////////////Clipping///////////////
			if (amplitude[OFDM_sym_idx] <= A) {
				s_clip_real[OFDM_sym_idx] = s_real[OFDM_sym_idx];
				s_clip_imag[OFDM_sym_idx] = s_imag[OFDM_sym_idx];
			}
			else {
				s_clip_real[OFDM_sym_idx] = A * cos(phase);
				s_clip_imag[OFDM_sym_idx] = A * sin(phase);
			}



			///////////////Calculating the deviation///////////////
			c_real[OFDM_sym_idx] = s_clip_real[OFDM_sym_idx] - s_real[OFDM_sym_idx];
			c_imag[OFDM_sym_idx] = s_clip_imag[OFDM_sym_idx] - s_imag[OFDM_sym_idx];

		}



		///////////////DFT///////////////
		DFT(c_real, c_imag, C_real, C_imag);



		///////////////ACE-Constraints///////////////
		for (OFDM_sym_idx = 0; OFDM_sym_idx < L * N; OFDM_sym_idx++) {

			if (OFDM_sym_idx < N) {

				///////////////constraints on real part///////////////
				if (S_real[OFDM_sym_idx] == -3.0 / sqrt(10) && C_real[OFDM_sym_idx] > 0) C_real[OFDM_sym_idx] = 0;
				if (S_real[OFDM_sym_idx] == -1.0 / sqrt(10)) C_real[OFDM_sym_idx] = 0;
				if (S_real[OFDM_sym_idx] == 1.0 / sqrt(10)) C_real[OFDM_sym_idx] = 0;
				if (S_real[OFDM_sym_idx] == 3.0 / sqrt(10) && C_real[OFDM_sym_idx] < 0) C_real[OFDM_sym_idx] = 0;

				///////////////constraints on imag part///////////////
				if (S_imag[OFDM_sym_idx] == -3.0 / sqrt(10) && C_imag[OFDM_sym_idx] > 0) C_imag[OFDM_sym_idx] = 0;
				if (S_imag[OFDM_sym_idx] == -1.0 / sqrt(10)) C_imag[OFDM_sym_idx] = 0;
				if (S_imag[OFDM_sym_idx] == 1.0 / sqrt(10)) C_imag[OFDM_sym_idx] = 0;
				if (S_imag[OFDM_sym_idx] == 3.0 / sqrt(10) && C_imag[OFDM_sym_idx] < 0) C_imag[OFDM_sym_idx] = 0;

			}
			else {
				C_real[OFDM_sym_idx] = 0;
				C_imag[OFDM_sym_idx] = 0;
			}

		}



		///////////////IDFT///////////////
		IDFT(C_real, C_imag, c_real, c_imag);



		///////////////Determining the Step Size///////////////
		for (OFDM_sym_idx = 0; OFDM_sym_idx < L * N; OFDM_sym_idx++) {
			c_proj[OFDM_sym_idx] = (s_real[OFDM_sym_idx] * c_real[OFDM_sym_idx] + s_imag[OFDM_sym_idx] * c_imag[OFDM_sym_idx]) / amplitude[OFDM_sym_idx];
		}

		Mu = 1000;
		for (OFDM_sym_idx = 0; OFDM_sym_idx < L * N; OFDM_sym_idx++) {
			delta[OFDM_sym_idx] = 0;
			if (c_proj[OFDM_sym_idx] > 0) {
				delta[OFDM_sym_idx] = (E - amplitude[OFDM_sym_idx]) / (c_proj[OFDM_sym_idx] - c_proj[n_max]);
				if (Mu > delta[OFDM_sym_idx]) Mu = delta[OFDM_sym_idx];
			}
		}





		if (Mu < 0 || Mu == 1000) {
			//cout << Mu; system("pause");
			Mu = 0;
		}
		///////////////New OFDM signal///////////////
		for (OFDM_sym_idx = 0; OFDM_sym_idx < L * N; OFDM_sym_idx++) {
			s_real[OFDM_sym_idx] = s_real[OFDM_sym_idx] + Mu * c_real[OFDM_sym_idx];
			s_imag[OFDM_sym_idx] = s_imag[OFDM_sym_idx] + Mu * c_imag[OFDM_sym_idx];
		}
	}





}


