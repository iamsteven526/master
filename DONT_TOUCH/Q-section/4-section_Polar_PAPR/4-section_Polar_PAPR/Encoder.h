
//using namespace std;


void Encoder(int Inf[k_c], int CodeWord[n_c], int G[k_c][n_c]) {

	int bit_idx, code_idx;



	for (code_idx = 0; code_idx < n_c; code_idx++) {
		CodeWord[code_idx] = 0;
	}

	for (bit_idx = 0; bit_idx < k_c; bit_idx++) {
		if (Inf[bit_idx] != 0) {
			for (code_idx = 0; code_idx < n_c; code_idx++) {
				CodeWord[code_idx] = CodeWord[code_idx] ^ G[bit_idx][code_idx];
			}
		}
	}


}







