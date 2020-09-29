#define		_CRT_SECURE_NO_WARNINGS


#include	<stdio.h>
#include	<fstream>
#include	<stdlib.h>
#include	<math.h>
#include	<random>
#include	<iostream>
#include	<vector>
#include	<iomanip>
//#include	<windows.h>



//using namespace std;

#define		Mod_file	 "16QAM_S2B_list.txt"
#define		log_name	 "LOG_Coded OFDM.txt"



#define		SNR_min					5
#define		SNR_max					5
#define		SNR_step				0.5


#define		PAPR_min				6
#define		PAPR_max				11
#define		PAPR_step				0.25
#define		tech_num				4
#define		S						16					// S=2^Q, Q-section, Q=4



#define		iter_max				60					// maximum number of LDPC decoding iterations
#define		iter_max_q				30					// maximum number of NBLDPC decoding iterations
#define		code_num_max			100000



#define		crc_num					0
#define		coderate				0.5


#define		n_c						1024				    // number of polar code bits
#define		k_c						512					    // number of polar information bits



#define		M						16					// Modulation order is equal to 16, 16-QAM
#define		p						4					// bits number of a modulation symbol, p = M^(1/2)


#define		N						256					// number of subcarriers
#define		L						4					// over-sampling factor


#define		PI						3.1415926			







#include	"AWGN.h"
#include	"DFT_IDFT.h"
#include	"Encoder.h"
#include	"LLR.h"
#include	"polar.h"
#include	"parameters.h"
#include	"ACE-SGP.h"





/*int			G_Polar[k_c][n_c] = { {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
							  {1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0},
							  {1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0},
							  {1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0} };*/


double		QAM_Constellation[M][2];			    //[0]: amplitude of real part, [1]: amplitude of imaginary part
int			QAM_S2B_Table[M][p];



int			row_idx, col_idx;
int			q_idx, s_idx;
double		SNR_idx;
double		PAPR_idx;
int			tech_idx, pattern_idx;
int			PAPR_num[tech_num];
int			sample_num = 1;
int			code_idx, bit_idx, sym_idx, mod_sym_idx, OFDM_sym_idx;
int			mod_idx;
double		snr, sigma_c[tech_num], sigmasquare[tech_num];




int			Inf_b[tech_num][k_c];
int			CodeWord[tech_num][n_c];
int			Mod[tech_num][N];
double		Mod_real[tech_num][N];
double		Mod_imag[tech_num][N];
double		temp_Ltimes_Mod_real[L * N];
double		temp_Ltimes_Mod_imag[L * N];
double		Ltimes_Mod_real[tech_num][L * N];
double		Ltimes_Mod_imag[tech_num][L * N];
double		temp_OFDM_real[L * N];
double		temp_OFDM_imag[L * N];
double		OFDM_real[tech_num][L * N];
double		OFDM_imag[tech_num][L * N];


double		amplitude[tech_num][L * N];
double		sita = 0;
double		phase = 0;
double		A_8dB;


double		temp_power;
double		avg_OFDM_power[2];			//[0]:origin, [1]:ACE-SGP
double		max_OFDM_power[2];			//[0]:origin, [1]:ACE-SGP
double		PAPR[tech_num];
double		min_OFDM_PAPR;
int			PAPR_sample[tech_num][5 + 4 * (PAPR_max - PAPR_min - 1)];
double		CCDF[tech_num][5 + 4 * (PAPR_max - PAPR_min - 1)];





double		noise[2];									// [0]: real part, [1]: imag part;
double		Received_OFDM_real[tech_num][L * N];
double		Received_OFDM_imag[tech_num][L * N];
double		Recovered_Mod_real[tech_num][L * N];
double		Recovered_Mod_imag[tech_num][L * N];
double		Priori_LLR[tech_num][N][p];
double		Polar_LLR[tech_num][n_c];
//std::vector<std::vector<int>> decoded_result(NUM_USER, std::vector<int>(DATA_LEN));
std::vector<int> decoded_result(k_c);
int			recovered_Inf_b[tech_num][k_c];
int			NBC_index[16];						//16-section, 16 bits
int			flag;




double		SER[tech_num], BER[tech_num], FER[tech_num];
int			Err_S_Num[tech_num], Err_B_Num[tech_num], Err_F_Num[tech_num], Error[tech_num];


int			cursor_idx;


int main() {
	

	srand(time(NULL));


	FILE* S2B_list;
	FILE* Log;






	///////////////QAM Constellation Construction///////////////
	S2B_list = fopen(Mod_file, "r");
	for (mod_sym_idx = 0; mod_sym_idx < M; mod_sym_idx++) {
		for (bit_idx = 0; bit_idx < p; bit_idx++) {
			fscanf(S2B_list, "%d", &QAM_S2B_Table[mod_sym_idx][bit_idx]);
		}
	}
	fclose(S2B_list);



	QAM_Constellation[0][0] = -3.0 / sqrt(10);								//average power normalized to 1
	QAM_Constellation[0][1] = 3.0 / sqrt(10);

	QAM_Constellation[1][0] = -3.0 / sqrt(10);
	QAM_Constellation[1][1] = 1.0 / sqrt(10);

	QAM_Constellation[2][0] = -3.0 / sqrt(10);
	QAM_Constellation[2][1] = -1.0 / sqrt(10);

	QAM_Constellation[3][0] = -3.0 / sqrt(10);
	QAM_Constellation[3][1] = -3.0 / sqrt(10);

	QAM_Constellation[4][0] = -1.0 / sqrt(10);
	QAM_Constellation[4][1] = -3.0 / sqrt(10);

	QAM_Constellation[5][0] = -1.0 / sqrt(10);
	QAM_Constellation[5][1] = -1.0 / sqrt(10);

	QAM_Constellation[6][0] = -1.0 / sqrt(10);
	QAM_Constellation[6][1] = 1.0 / sqrt(10);

	QAM_Constellation[7][0] = -1.0 / sqrt(10);
	QAM_Constellation[7][1] = 3.0 / sqrt(10);

	QAM_Constellation[8][0] = 1.0 / sqrt(10);
	QAM_Constellation[8][1] = 3.0 / sqrt(10);

	QAM_Constellation[9][0] = 1.0 / sqrt(10);
	QAM_Constellation[9][1] = 1.0 / sqrt(10);

	QAM_Constellation[10][0] = 1.0 / sqrt(10);
	QAM_Constellation[10][1] = -1.0 / sqrt(10);

	QAM_Constellation[11][0] = 1.0 / sqrt(10);
	QAM_Constellation[11][1] = -3.0 / sqrt(10);

	QAM_Constellation[12][0] = 3.0 / sqrt(10);
	QAM_Constellation[12][1] = -3.0 / sqrt(10);

	QAM_Constellation[13][0] = 3.0 / sqrt(10);
	QAM_Constellation[13][1] = -1.0 / sqrt(10);

	QAM_Constellation[14][0] = 3.0 / sqrt(10);
	QAM_Constellation[14][1] = 1.0 / sqrt(10);

	QAM_Constellation[15][0] = 3.0 / sqrt(10);
	QAM_Constellation[15][1] = 3.0 / sqrt(10);













	///////////////START///////////////
	Log = fopen(log_name, "a");
	fprintf(Log, "iter_max=%d\n", iter_max);
	fclose(Log);




	for (SNR_idx = SNR_min; SNR_idx <= SNR_max; SNR_idx += SNR_step) {

		PolarCode polar(10, k_c, crc_num, 1);
		polar.initialize_frozen_bits(0, 1);
		snr = pow(10, 0.1 * SNR_idx);
		sigma_c[0] = sqrt(1 / (2 * (double(k_c + crc_num) / double(n_c)) * snr * p));
		sigmasquare[0] = sigma_c[0] * sigma_c[0];

		sigma_c[1] = sigma_c[0]; sigmasquare[1] = sigmasquare[0];

		sigma_c[2] = sqrt(1 / (2 * (double(k_c + crc_num) / double(n_c)) * snr * p));
		sigmasquare[2] = sigma_c[2] * sigma_c[2];

		sigma_c[3] = sigma_c[2]; sigmasquare[3] = sigmasquare[2];


		for (tech_idx = 0; tech_idx < tech_num; tech_idx++) {
			Err_B_Num[tech_idx] = 0;
			Err_F_Num[tech_idx] = 0;
			Err_S_Num[tech_idx] = 0;
			Error[tech_idx] = 0;
		}
		for (code_idx = 0; code_idx < code_num_max; code_idx++) {




			///////////////Bit-stream Generating///////////////
			for (bit_idx = 0; bit_idx < k_c; bit_idx++) {

				if (rand() > RAND_MAX / 2) Inf_b[0][bit_idx] = 1;
				else                       Inf_b[0][bit_idx] = 0;
			}
			//for (bit_idx = 0; bit_idx < k_c; bit_idx++) std::cout << Inf_b[bit_idx] << " "; std::cout << std::endl << std::endl; system("pause");
			//for (bit_idx = 0; bit_idx < k_c; bit_idx++) Inf_b[1][bit_idx] = Inf_b[0][bit_idx];



			///////////////find the NBC index///////////////
			/*int test[n_c];
			for (bit_idx = 0; bit_idx < n_c; bit_idx++) test[bit_idx] = 0;


			///////////////v0///////////////
			//for (bit_idx = 0; bit_idx < n_c; bit_idx++) test[bit_idx] = 1;


			///////////////v1///////////////
			//for (bit_idx = 0; bit_idx < 512; bit_idx++) test[bit_idx] = 1;


			///////////////v2///////////////
			//for (bit_idx = 0; bit_idx < 256; bit_idx++) test[bit_idx] = 1;
			//for (bit_idx = 512; bit_idx < 768; bit_idx++) test[bit_idx] = 1;


			///////////////v3///////////////
			//for (bit_idx = 0; bit_idx < 128; bit_idx++) test[bit_idx] = 1;
			//for (bit_idx = 256; bit_idx < 384; bit_idx++) test[bit_idx] = 1;
			//for (bit_idx = 512; bit_idx < 640; bit_idx++) test[bit_idx] = 1;
			//for (bit_idx = 768; bit_idx < 896; bit_idx++) test[bit_idx] = 1;


			///////////////v4///////////////
			//for (bit_idx = 0; bit_idx < 64; bit_idx++) test[bit_idx] = 1;
			//for (bit_idx = 128; bit_idx < 192; bit_idx++) test[bit_idx] = 1;
			//for (bit_idx = 256; bit_idx < 320; bit_idx++) test[bit_idx] = 1;
			//for (bit_idx = 384; bit_idx < 448; bit_idx++) test[bit_idx] = 1;
			//for (bit_idx = 512; bit_idx < 576; bit_idx++) test[bit_idx] = 1;
			//for (bit_idx = 640; bit_idx < 704; bit_idx++) test[bit_idx] = 1;
			//for (bit_idx = 768; bit_idx < 832; bit_idx++) test[bit_idx] = 1;
			//for (bit_idx = 896; bit_idx < 960; bit_idx++) test[bit_idx] = 1;


			int show[n_c];
			///////////////Polar Encoding///////////////
			for (row_idx = 0; row_idx < k_c; row_idx++) {
				//std::cout << row_idx << std::endl;
				for (bit_idx = 0; bit_idx < k_c; bit_idx++) Inf_b[0][bit_idx] = 0;
				Inf_b[0][row_idx] = 1;
				polar.encode(Inf_b[0], CodeWord[0], 0);
				for (bit_idx = 0; bit_idx < n_c; bit_idx++) show[bit_idx] = (test[bit_idx] + CodeWord[0][bit_idx]) % 2;
				//std::cout << row_idx << ". ";
				//for (bit_idx = 0; bit_idx < n_c; bit_idx++) std::cout << show[bit_idx] << " ";
				//std::cout << std::endl << std::endl;
				for (bit_idx = 0; bit_idx < n_c; bit_idx++) {
					if (show[bit_idx] != 0) break;
				}
				if (bit_idx == n_c) {
					std::cout << row_idx;
					system("pause");
				}
			}
			system("pause");*/

			NBC_index[0] = 0; NBC_index[1] = 1; NBC_index[2] = 2; NBC_index[3] = 3;
			NBC_index[4] = 4; NBC_index[5] = 11; NBC_index[6] = 12; NBC_index[7] = 14;
			NBC_index[8] = 13; NBC_index[9] = 15; NBC_index[10] = 16; NBC_index[11] = 53;
			NBC_index[12] = 54; NBC_index[13] = 55; NBC_index[14] = 56; NBC_index[15] = 147;




			//for (bit_idx = 0; bit_idx < k_c; bit_idx++) Inf_b[0][bit_idx] = 0; Inf_b[0][147] = 1;
			polar.encode(Inf_b[0], CodeWord[0], 0);
			//for (bit_idx = 0; bit_idx < n_c; bit_idx++) std::cout << CodeWord[0][bit_idx] << " "; system("pause");



			///////////////16-QAM Modulation with Gray Mapping///////////////
			for (mod_sym_idx = 0; mod_sym_idx < N; mod_sym_idx++) {


				for (mod_idx = 0; mod_idx < M; mod_idx++) {
					for (bit_idx = 0; bit_idx < p; bit_idx++) {
						if (CodeWord[0][bit_idx + mod_sym_idx * p] != QAM_S2B_Table[mod_idx][bit_idx]) break;
					}
					if (bit_idx == p) {
						Mod[0][mod_sym_idx] = mod_idx;
						Mod_real[0][mod_sym_idx] = QAM_Constellation[mod_idx][0];
						Mod_imag[0][mod_sym_idx] = QAM_Constellation[mod_idx][1];
						break;
					}
				}

			}
			/*for (mod_sym_idx = 0; mod_sym_idx < N; mod_sym_idx++) cout << mod_sym_idx << ".  " << Mod[mod_sym_idx] << endl;
			cout << endl << endl << endl; system("pause");*/




			///////////////Over-Sampling///////////////
			for (mod_sym_idx = 0; mod_sym_idx < L * N; mod_sym_idx++) {
				if (mod_sym_idx < N) {
					Ltimes_Mod_real[0][mod_sym_idx] = Mod_real[0][mod_sym_idx];
					Ltimes_Mod_imag[0][mod_sym_idx] = Mod_imag[0][mod_sym_idx];
				}
				else {
					Ltimes_Mod_real[0][mod_sym_idx] = 0;
					Ltimes_Mod_imag[0][mod_sym_idx] = 0;
				}
			}





			///////////////OFDM Signal Modulation///////////////
			IDFT(Ltimes_Mod_real[0], Ltimes_Mod_imag[0], OFDM_real[0], OFDM_imag[0]);

			//for (bit_idx = 0; bit_idx < n_c; bit_idx++) std::cout << OFDM_real[0][bit_idx] << " "; system("pause");
			//std::cout << std::endl << std::endl;


			///////////////ACE-SGP PAPR Reduction///////////////
			ACE_SGP(Ltimes_Mod_real[0], Ltimes_Mod_imag[0], OFDM_real[1], OFDM_imag[1]);

			//for (bit_idx = 0; bit_idx < n_c; bit_idx++) std::cout << OFDM_real[1][bit_idx] << " "; system("pause");
			//std::cout << std::endl << std::endl;



			///////////////PAPR calculation///////////////
			for (tech_idx = 0; tech_idx < 2; tech_idx++) {
				max_OFDM_power[tech_idx] = 0;
				avg_OFDM_power[tech_idx] = 0;
			}
			for (OFDM_sym_idx = 0; OFDM_sym_idx < L * N; OFDM_sym_idx++) {
				for (tech_idx = 0; tech_idx < 2; tech_idx++) {
					temp_power = pow(OFDM_real[tech_idx][OFDM_sym_idx], 2) + pow(OFDM_imag[tech_idx][OFDM_sym_idx], 2);
					avg_OFDM_power[tech_idx] += temp_power;
					if (temp_power > max_OFDM_power[tech_idx]) {
						max_OFDM_power[tech_idx] = temp_power;
					}
				}
			}
			for (tech_idx = 0; tech_idx < 2; tech_idx++) avg_OFDM_power[tech_idx] = avg_OFDM_power[tech_idx] / (L * N);
			/*if (max_OFDM_s_power[0] != max_OFDM_s_power[1]) {
				cout <<  max_OFDM_s_power[0] << " " << max_OFDM_s_power[1];
				system("pause");
			}*/
			for (tech_idx = 0; tech_idx < 2; tech_idx++) PAPR[tech_idx] = 10 * log10(max_OFDM_power[tech_idx] / avg_OFDM_power[tech_idx]);

			//std::cout << PAPR[0] << " " << PAPR[1]; system("pause");



			///////////////Q-Section PAPR Reduction///////////////
			int digits[4] = { 0 }, index = 0;
			int quotient, remainder, xx;
			min_OFDM_PAPR = 100;
			for (int x = 0; x < 16; x++) {
				for (bit_idx = 0; bit_idx < k_c; bit_idx++) Inf_b[1][bit_idx] = Inf_b[0][bit_idx];
				for (bit_idx = 0; bit_idx < 4; bit_idx++) digits[bit_idx] = 0;
				index = 0;
				xx = x;
				while (xx != 0) {
					quotient = xx / 2;
					remainder = xx % 2;
					if (remainder == 1) digits[index++] = 1;
					else digits[index++] = 0;
					xx = xx / 2;
				}
				if (digits[0] == 1) Inf_b[1][0] = (Inf_b[1][0] + 1) % 2;
				if (digits[1] == 1) Inf_b[1][1] = (Inf_b[1][1] + 1) % 2;
				if (digits[2] == 1) Inf_b[1][2] = (Inf_b[1][2] + 1) % 2;
				if (digits[3] == 1) Inf_b[1][11] = (Inf_b[1][11] + 1) % 2;
				///////////////Polar Encoding///////////////
				polar.encode(Inf_b[1], CodeWord[1], 0);
				///////////////16-QAM Modulation with Gray Mapping///////////////
				for (mod_sym_idx = 0; mod_sym_idx < N; mod_sym_idx++) {


					for (mod_idx = 0; mod_idx < M; mod_idx++) {
						for (bit_idx = 0; bit_idx < p; bit_idx++) {
							if (CodeWord[1][bit_idx + mod_sym_idx * p] != QAM_S2B_Table[mod_idx][bit_idx]) break;
						}
						if (bit_idx == p) {
							Mod[1][mod_sym_idx] = mod_idx;
							Mod_real[1][mod_sym_idx] = QAM_Constellation[mod_idx][0];
							Mod_imag[1][mod_sym_idx] = QAM_Constellation[mod_idx][1];
							break;
						}
					}

				}
				///////////////Over-Sampling///////////////
				for (mod_sym_idx = 0; mod_sym_idx < L * N; mod_sym_idx++) {
					if (mod_sym_idx < N) {
						temp_Ltimes_Mod_real[mod_sym_idx] = Mod_real[1][mod_sym_idx];
						temp_Ltimes_Mod_imag[mod_sym_idx] = Mod_imag[1][mod_sym_idx];
					}
					else {
						temp_Ltimes_Mod_real[mod_sym_idx] = 0;
						temp_Ltimes_Mod_imag[mod_sym_idx] = 0;
					}
				}
				///////////////OFDM Signal Modulation///////////////
				IDFT(temp_Ltimes_Mod_real, temp_Ltimes_Mod_imag, temp_OFDM_real, temp_OFDM_imag);
				///////////////PAPR calculation///////////////
				max_OFDM_power[0] = 0;	avg_OFDM_power[0] = 0;
				for (OFDM_sym_idx = 0; OFDM_sym_idx < L * N; OFDM_sym_idx++) {
					temp_power = pow(temp_OFDM_real[OFDM_sym_idx], 2) + pow(temp_OFDM_imag[OFDM_sym_idx], 2);
					avg_OFDM_power[0] += temp_power;
					if (temp_power > max_OFDM_power[0]) {
						max_OFDM_power[0] = temp_power;
					}
				}
				avg_OFDM_power[0] = avg_OFDM_power[0] / (L * N);
				/*if (max_OFDM_s_power[0] != max_OFDM_s_power[1]) {
					cout <<  max_OFDM_s_power[0] << " " << max_OFDM_s_power[1];
					system("pause");
				}*/
				if (min_OFDM_PAPR > 10 * log10(max_OFDM_power[0] / avg_OFDM_power[0])) {
					for (OFDM_sym_idx = 0; OFDM_sym_idx < L * N; OFDM_sym_idx++) {
						Ltimes_Mod_real[1][OFDM_sym_idx] = temp_Ltimes_Mod_real[OFDM_sym_idx];
						Ltimes_Mod_imag[1][OFDM_sym_idx] = temp_Ltimes_Mod_imag[OFDM_sym_idx];
						OFDM_real[2][OFDM_sym_idx] = temp_OFDM_real[OFDM_sym_idx];
						OFDM_imag[2][OFDM_sym_idx] = temp_OFDM_imag[OFDM_sym_idx];
					}
					min_OFDM_PAPR = 10 * log10(max_OFDM_power[0] / avg_OFDM_power[0]);
				}
			}
			PAPR[2] = min_OFDM_PAPR;
			//std::cout << PAPR[1] << std::endl; system("pause");


			///////////////ACE-SGP PAPR Reduction///////////////
			ACE_SGP(Ltimes_Mod_real[1], Ltimes_Mod_imag[1], OFDM_real[3], OFDM_imag[3]);
			///////////////PAPR calculation///////////////
			max_OFDM_power[0] = 0;	avg_OFDM_power[0] = 0;
			for (OFDM_sym_idx = 0; OFDM_sym_idx < L * N; OFDM_sym_idx++) {
				temp_power = pow(OFDM_real[3][OFDM_sym_idx], 2) + pow(OFDM_imag[3][OFDM_sym_idx], 2);
				avg_OFDM_power[0] += temp_power;
				if (temp_power > max_OFDM_power[0]) {
					max_OFDM_power[0] = temp_power;
				}
			}
			avg_OFDM_power[0] = avg_OFDM_power[0] / (L * N);
			PAPR[3] = 10 * log10(max_OFDM_power[0] / avg_OFDM_power[0]);




			///////////////Non-linear channel of PA///////////////
			for (tech_idx = 0; tech_idx < tech_num; tech_idx++) {
				for (OFDM_sym_idx = 0; OFDM_sym_idx < L * N; OFDM_sym_idx++) {
					amplitude[tech_idx][OFDM_sym_idx] = sqrt(pow(OFDM_real[tech_idx][OFDM_sym_idx], 2) + pow(OFDM_imag[tech_idx][OFDM_sym_idx], 2));
					sita = atan(OFDM_real[tech_idx][OFDM_sym_idx] / OFDM_imag[tech_idx][OFDM_sym_idx]);
					if (OFDM_real[tech_idx][OFDM_sym_idx] < 0) phase = sita + PI;
					else phase = sita;

					///////////////Clipping///////////////

					//A_8dB is the amplititude of saturation region for OBO = 8 dB

					A_8dB = sqrt(pow(10, 0.8)*(1 / double(L)));
					if (amplitude[tech_idx][OFDM_sym_idx] > A_8dB) {
						OFDM_real[tech_idx][OFDM_sym_idx] = A_8dB * cos(phase);
						OFDM_imag[tech_idx][OFDM_sym_idx] = A_8dB * sin(phase);
					}
					else {
						OFDM_real[tech_idx][OFDM_sym_idx] = OFDM_real[tech_idx][OFDM_sym_idx];
						OFDM_imag[tech_idx][OFDM_sym_idx] = OFDM_imag[tech_idx][OFDM_sym_idx];
					}

				}
			}
			




			///////////////AWGN///////////////
			for (OFDM_sym_idx = 0; OFDM_sym_idx < L * N; OFDM_sym_idx++) {
				noise[0] = normal(); noise[1] = normal();
				for (tech_idx = 0; tech_idx < tech_num; tech_idx++) {
					Received_OFDM_real[tech_idx][OFDM_sym_idx] = OFDM_real[tech_idx][OFDM_sym_idx] + noise[0] * sigma_c[tech_idx];
					Received_OFDM_imag[tech_idx][OFDM_sym_idx] = OFDM_imag[tech_idx][OFDM_sym_idx] + noise[1] * sigma_c[tech_idx];
				}
			}





			///////////////OFDM Signal DeModulation///////////////
			for (tech_idx = 0; tech_idx < tech_num; tech_idx++) {
				DFT(Received_OFDM_real[tech_idx], Received_OFDM_imag[tech_idx], Recovered_Mod_real[tech_idx], Recovered_Mod_imag[tech_idx]);
			}







			///////////////16-QAM DeModulation///////////////
			for (tech_idx = 0; tech_idx < tech_num; tech_idx++) {
				for (mod_sym_idx = 0; mod_sym_idx < N; mod_sym_idx++) {
					LLR(Recovered_Mod_real[tech_idx][mod_sym_idx], Recovered_Mod_imag[tech_idx][mod_sym_idx], sigmasquare[tech_idx], Priori_LLR[tech_idx][mod_sym_idx]);
				}
				for (bit_idx = 0; bit_idx < n_c; bit_idx++) Polar_LLR[tech_idx][bit_idx] = Priori_LLR[tech_idx][bit_idx / p][bit_idx % p];
			}
			/*for (tech_idx = 0; tech_idx < tech_num; tech_idx++) {
				for (bit_idx = 0; bit_idx < n_c; bit_idx++) {
					std::cout << Polar_LLR[tech_idx][bit_idx] << " ";
				}
				std::cout << std::endl << std::endl;
			}
			system("pause");*/




			///////////////Polar decoding///////////////
			for (tech_idx = 0; tech_idx < tech_num; tech_idx++) {
				decoded_result = polar.decode_scl_llr(LIST_SIZE, Polar_LLR[tech_idx], 0);
				for (bit_idx = 0; bit_idx < k_c; bit_idx++) {
					recovered_Inf_b[tech_idx][bit_idx] = decoded_result[bit_idx];
				}
			}

			//for (bit_idx = 0; bit_idx < k_c; bit_idx++) std::cout << Recovered_Inf_b[0][bit_idx] << " "; std::cout << std::endl << std::endl; system("pause");



            /*
			HANDLE hOut;
			COORD OutChar;
			OutChar.X = 0;
			OutChar.Y = cursor_idx * 6;
			hOut = GetStdHandle(STD_OUTPUT_HANDLE);
			SetConsoleCursorPosition(hOut, OutChar);
			*/
			///////////////Error Performance Measurements///////////////
			for (bit_idx = 0; bit_idx < k_c; bit_idx++) {
				flag = 0;
				for (int NBC_idx = 0; NBC_idx < 16; NBC_idx++) {
					if (bit_idx == NBC_index[NBC_idx]) {
						flag = 1;
						break;
					}
				}
				if (flag == 0) {
					for (tech_idx = 0; tech_idx < tech_num; tech_idx++) {
						if (Inf_b[0][bit_idx] != recovered_Inf_b[tech_idx][bit_idx]) {
							Error[tech_idx] = 1;
							Err_B_Num[tech_idx]++;
						}
					}
				}
			}
			for (tech_idx = 0; tech_idx < tech_num; tech_idx++) {
				Err_F_Num[tech_idx] += Error[tech_idx];
				Error[tech_idx] = 0;
			}


			std::cout << "SNR: " << std::setw(2) << SNR_idx;
			std::cout << ", frame = " << code_idx + 1 << std::flush;
			std::cout << "\n";
			for (tech_idx = 0; tech_idx < tech_num; tech_idx++) {
				BER[tech_idx] = (double)Err_B_Num[tech_idx] / double(code_idx + 1) / (k_c - 16);
				FER[tech_idx] = (double)Err_F_Num[tech_idx] / double(code_idx + 1);
				std::cout << "Error_frame = " << Err_F_Num[tech_idx] << std::flush;
				std::cout << ", BER  : " << BER[tech_idx] << std::flush;
				std::cout << ", FER  : " << FER[tech_idx] << std::flush;
				std::cout << "\n";
			}





            /*
			OutChar.X = 0;
			OutChar.Y = 80;
			hOut = GetStdHandle(STD_OUTPUT_HANDLE);
			SetConsoleCursorPosition(hOut, OutChar);
			*/
			///////////////PAPR Performance Measurements///////////////
			for (tech_idx = 0; tech_idx < tech_num; tech_idx++) PAPR_num[tech_idx] = 0;
			for (PAPR_idx = PAPR_min; PAPR_idx <= PAPR_max; PAPR_idx += PAPR_step) {
				for (tech_idx = 0; tech_idx < tech_num; tech_idx++) {
					if (PAPR[tech_idx] > PAPR_idx) PAPR_sample[tech_idx][PAPR_num[tech_idx]]++;
					CCDF[tech_idx][PAPR_num[tech_idx]] = double(PAPR_sample[tech_idx][PAPR_num[tech_idx]]) / sample_num;
					PAPR_num[tech_idx]++;
				}
			}



			PAPR_num[0] = 0;
			for (PAPR_idx = PAPR_min; PAPR_idx <= PAPR_max; PAPR_idx += PAPR_step) {
				std::cout << "PAPR(dB): " << PAPR_idx << ", " << "no PAPR reduction_" << CCDF[0][PAPR_num[0]] << "   " << "ACE-SGP_" << CCDF[1][PAPR_num[0]] << "   " << "Q-PTS_" << CCDF[2][PAPR_num[0]] << "   " << "Q-PTS-ACE-SGP_" << CCDF[3][PAPR_num[0]] << std::flush;
				std::cout << "\n";
				PAPR_num[0]++;
			}

			sample_num++;
		}

		cursor_idx++;
	}




	system("pause");

	return 0;
}
