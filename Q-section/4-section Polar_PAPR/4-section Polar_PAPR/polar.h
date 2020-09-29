

#include <iostream>
#include <vector>
#include <stack>

class PolarCode
{

public:

	PolarCode(uint8_t num_layers, uint16_t info_length, uint16_t crc_size, int num_users)
	{
		_n = num_layers;
		_info_length = info_length;
		_crc_size = crc_size;
		_llr_based_computation = true;
		_block_length = (uint16_t)(1 << _n);
		_crc_matrix.resize(_crc_size);
		_bit_rev_order.resize(_block_length);
		create_bit_rev_order();

		_frozen_bits.resize(num_users);
		_channel_order_sorted.resize(num_users);
		_channel_vec.resize(num_users);
		for (int nuser = 0; nuser < num_users; nuser++)
		{
			_frozen_bits.at(nuser).resize(_block_length);
			_channel_order_sorted.at(nuser).resize(_block_length);
			_channel_vec.at(nuser).resize(_block_length);
		}

		G.resize(_block_length);
		for (int row = 0; row < _block_length; row++)
		{
			G.at(row).resize(_block_length);
		}
	}

	void initialize_frozen_bits(uint8_t type, double designedSNR);
	void encode(int *info_bits, int *codeword, int nuser);
	void  encode_bp(int* info_bits, int* codeword, int nuser);
	std::vector<int> decode_scl_llr(uint16_t list_size, double *llr, int nuser);
	void decode_jpscl_llr(uint16_t list_size, double **app, std::vector<std::vector<int>>& decodedResult); // JPSCL
	void decode_BP_Joint(double **app, std::vector<std::vector<int>>& detectedResult);
	std::vector<int> decode_BP_sep(double *llr, int nuser);
	std::vector<double> vector_product(std::vector<double> v1, std::vector<double> v2);

private:

	uint8_t _n;
	uint16_t _info_length;
	uint16_t _block_length;
	uint16_t _crc_size;

	double _design_parameter;

	std::vector<std::vector<double>> _channel_vec;
	std::vector<std::vector<uint8_t>> _frozen_bits;
	std::vector<std::vector<uint16_t>> _channel_order_sorted;
	std::vector<std::vector<uint8_t>> _crc_matrix;
	std::vector<uint16_t> _bit_rev_order;
	std::vector<std::vector<int>> G;

	std::vector<uint16_t> puncturing_pattern;

	void create_bit_rev_order();

	double JFunction(double sigma);
	double InverseJFunction(double I);
	double min_sum(double a, double b);

	std::vector<uint8_t> decode_scl_p1(std::vector<double> p1, std::vector<double> p0, uint16_t list_size, int nuser);
	std::vector<int> decode_scl(int nuser);
	bool _llr_based_computation;

	std::vector<std::vector<double *>> _arrayPointer_LLR;
	std::vector<std::vector<double**>> _arrayPointer_pr_vector; // JPSCL
	std::vector<double> _pathMetric_LLR;

	uint16_t _list_size;

	std::stack<uint16_t> _inactivePathIndices;
	std::vector<uint16_t > _activePath;

	std::vector<std::vector<double*>> _arrayPointer_P;
	std::vector<std::vector<uint8_t*>> _arrayPointer_C;
	std::vector<std::vector<double**>> _arrayPointer_C_vector; // JPSCL

	std::vector<uint8_t *> _arrayPointer_Info;
	std::vector<std::vector<uint16_t>> _pathIndexToArrayIndex;
	std::vector<std::stack<uint16_t>> _inactiveArrayIndices;
	std::vector<std::vector<uint16_t>> _arrayReferenceCount;


	void initializeDataStructures();
	uint16_t assignInitialPath();
	uint16_t clonePath(uint16_t);
	void killPath(uint16_t l);

	double *getArrayPointer_P(uint16_t lambda, uint16_t  l);
	double *getArrayPointer_LLR(uint16_t lambda, uint16_t  l);
	double **getArrayPointer_pr_vector(uint16_t lambda, uint16_t  l);
	uint8_t *getArrayPointer_C(uint16_t lambda, uint16_t  l);
	double **getArrayPointer_C_vector(uint16_t lambda, uint16_t  l);

	void recursivelyCalcP(uint16_t lambda, uint16_t phi);
	void recursivelyCalcLLR(uint16_t lambda, uint16_t phi);
	void recursivelyUpdateC(uint16_t lambda, uint16_t phi);

	void continuePaths_FrozenBit(uint16_t phi, int nuser);
	void continuePaths_UnfrozenBit(uint16_t phi, int nuser);

	uint16_t findMostProbablePath(bool check_crc, int nuser);

	bool crc_check(uint8_t * info_bits_padded, int nuser);

};


