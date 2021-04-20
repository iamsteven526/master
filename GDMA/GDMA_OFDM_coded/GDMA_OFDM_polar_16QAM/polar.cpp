#include <iostream>
#include <random>
#include <algorithm>
#include <functional>
#include <fstream>
#include <algorithm>
#include "polar.h"
#include "parameters.h"
#include <cstring>
#include <math.h>

void PolarCode::initialize_frozen_bits(uint8_t type, double designedSNR)
{
	for (int nuser = 0; nuser < NUM_USER; nuser++)
	{
		if (type != 4 || type != 5)
		{
			if (type == 0) // Bhattacharyya bound
			{
				_design_parameter = exp(-pow(10., designedSNR / 10.));
			}
			else if (type == 1 || type == 2 || type == 3) // BI-AWGN capacity
			{
				double variance = 0.5 * pow(10., -designedSNR / 10.);
				_design_parameter = JFunction(2. / sqrt(variance));
			}
			for (uint16_t i = 0; i < _block_length; ++i)
			{
				_channel_vec.at(nuser).at(i) = _design_parameter;
			}
			for (uint8_t iteration = 0; iteration < _n; ++iteration)
			{
				uint16_t increment = 1 << iteration;
				for (uint16_t j = 0; j < increment; j += 1)
				{
					for (uint16_t i = 0; i < _block_length; i += 2 * increment)
					{
						double c1 = _channel_vec.at(nuser).at(i + j);
						double c2 = _channel_vec.at(nuser).at(i + j + increment);
						if (type == 0) // Bhattacharyya bound
						{
							_channel_vec.at(nuser).at(i + j) = c1 + c2 - c1 * c2;
							_channel_vec.at(nuser).at(i + j + increment) = c1 * c2;
						}
						if (type == 1) // Capacity bound
						{
							_channel_vec.at(nuser).at(i + j) = c1 * c2;
							_channel_vec.at(nuser).at(i + j + increment) = c1 + c2 - c1 * c2;
						}
						if (type == 2) // Piecewise linear approximation
						{
							if (c1 >= 0.1618)
							{
								_channel_vec.at(nuser).at(i + j + increment) = c1 + c1 * (-0.951 * c1 + 0.951);
							}
							else
							{
								_channel_vec.at(nuser).at(i + j + increment) = c1 + c1 * (-1.198 * c1 + 0.996);
							}
							_channel_vec.at(nuser).at(i + j) = 2 * c1 - _channel_vec.at(nuser).at(i + j + increment);
						}
						if (type == 3) // Gaussain approximation
						{
							_channel_vec.at(nuser).at(i + j) = 1. - JFunction(sqrt(pow(InverseJFunction(1 - c1), 2.) * 2.));
							_channel_vec.at(nuser).at(i + j + increment) = JFunction(sqrt(pow(InverseJFunction(c1), 2) * 2));
						}
					}
				}
			}
		}
		std::size_t n_t(0);
		std::generate(std::begin(_channel_order_sorted.at(nuser)), std::end(_channel_order_sorted.at(nuser)), [&] { return n_t++; });
		if (type == 0) // ascending order
		{
			std::sort(std::begin(_channel_order_sorted.at(nuser)),
				std::end(_channel_order_sorted.at(nuser)),
				[&](int i1, int i2) { return _channel_vec.at(nuser).at(_bit_rev_order.at(i1)) < _channel_vec.at(nuser).at(_bit_rev_order.at(i2)); });
		}
		else if (type == 1 || type == 2 || type == 3) // descending order
		{
			std::sort(std::begin(_channel_order_sorted.at(nuser)),
				std::end(_channel_order_sorted.at(nuser)),
				[&](int i1, int i2) { return _channel_vec.at(nuser).at(_bit_rev_order.at(i1)) > _channel_vec.at(nuser).at(_bit_rev_order.at(i2)); });
		}
		else if (type == 4)
		{
			std::fstream sequence;
			sequence.open("3gpp_polar_seq.txt", std::ios::in);
			for (int i = 0; i < 1024; i++) // fixed code length of 1024 from 3GPP definition 
			{
				sequence >> _channel_order_sorted.at(nuser).at(i);
				_channel_order_sorted.at(nuser).at(i) = 1023 - _channel_order_sorted.at(nuser).at(i);
			}
			sequence.close();
		}
		else if (type == 5)
		{
			std::fstream sequence;
			sequence.open("NBC_polar_seq.txt", std::ios::in);
			for (int i = 0; i < 1024; i++) // fixed code length of 1024 from 3GPP definition 
			{
				sequence >> _channel_order_sorted.at(nuser).at(i);
				_channel_order_sorted.at(nuser).at(i) = 1023 - _channel_order_sorted.at(nuser).at(i);
			}
			sequence.close();
		}
		else
		{
			printf("\nPARAMETER SETTING IS WRONG\n");
			//system("pause");
		}
		uint16_t effective_info_length = _info_length + _crc_size;
		for (uint16_t i = 0; i < effective_info_length; ++i)
		{
			_frozen_bits.at(nuser).at(_channel_order_sorted.at(nuser).at(i)) = 0;
		}
		for (uint16_t i = effective_info_length; i < _block_length; ++i)
		{
			_frozen_bits.at(nuser).at(_channel_order_sorted.at(nuser).at((i))) = 1;
		}
	}
	for (uint8_t bit = 0; bit < _crc_size; ++bit)
	{
		_crc_matrix.at(bit).resize(_info_length);
		for (uint16_t info_bit = 0; info_bit < _info_length; ++info_bit)
		{
			_crc_matrix.at(bit).at(info_bit) = (uint8_t)(rand() % 2);
		}
	}

	//creat G for encode_bp

	int n = log(_block_length) / log(2);
	int x;

	G[0][0] = 1;

	for (int i = 0; i < n; i++)
	{
		x = pow(2, i);
		for (int j = 0; j < pow(2, i); j++)
		{
			for (int k = 0; k < pow(2, i); k++)
			{
				G[j + x][k] = G[j][k];
				G[j + x][k + x] = G[j][k];
			}
		}
	}
	
	/*for (int i = 0; i < _block_length; i++)
	{
		for (int j = 0; j < _block_length; j++)
		{
			std::cout << G[i][j] << " ";
		}
		std::cout << std::endl;
	}*/
}

void PolarCode::encode(int *info_bits, int* codeword, int nuser)
{
	std::vector<uint8_t> info_bits_padded(_block_length, 0);
	for (uint16_t i = 0; i < _info_length; ++i)
	{
		info_bits_padded.at(_channel_order_sorted.at(nuser).at(i)) = info_bits[i];
	}
	for (uint16_t i = _info_length; i < _info_length + _crc_size; ++i)
	{
		uint8_t crc_bit = 0;
		for (uint16_t j = 0; j < _info_length; ++j)
		{
			crc_bit = (uint8_t)((crc_bit + _crc_matrix.at(i - _info_length).at(j) * info_bits[j]) % 2);
		}
		info_bits_padded.at(_channel_order_sorted.at(nuser).at(i)) = crc_bit;
	}
	for (uint8_t iteration = 0; iteration < _n; ++iteration)
	{
		uint16_t  increment = (uint16_t)(1 << iteration);
		for (uint16_t j = 0; j < increment; j += 1)
		{
			for (uint16_t i = 0; i < _block_length; i += 2 * increment)
			{
				info_bits_padded.at(i + j) = (uint8_t)((info_bits_padded.at(i + j) + info_bits_padded.at(i + j + increment)) % 2);
			}
		}
	}


	for (uint16_t i = 0; i < _block_length; ++i)
	{
		codeword[i] = info_bits_padded.at(_bit_rev_order.at(i));
	}

}

void PolarCode::encode_bp(int *info_bits, int* codeword, int nuser)
{
	sort(_channel_order_sorted.at(nuser).begin(), _channel_order_sorted.at(nuser).begin()+DATA_LEN);

	int reg;
	for (int i = 0; i < _block_length; i++)
	{
		reg = 0;
		for (int j = 0; j < DATA_LEN; j++)
		{
			reg = reg + G[_channel_order_sorted[nuser][j]][i] * info_bits[j];
			//	cout <<G[i][j];
			//cout << random_message[j] << " ";
		}
		//std::cout << std::endl;
		codeword[i] = reg % 2;
	}
}

bool PolarCode::crc_check(uint8_t *info_bit_padded, int nuser)
{
	std :: vector<std :: vector<uint8_t>> info_bit_binary((int8_t)NUM_USER ,std::vector<uint8_t>(_info_length + _crc_size));
	//std::vector<uint8_t> info_bit_binary(_info_length);
	for (uint16_t j = 0; j < _info_length + _crc_size; ++j)
	{
		uint8_t reg = 0;

		if (!JOINT)
		{
			info_bit_binary.at(0).at(j) = info_bit_padded[_channel_order_sorted.at(nuser).at(j)];
			for (int n = 1; n < NUM_USER; n++)
			{
				info_bit_binary.at(n).at(j) = info_bit_binary.at(0).at(j);
			}
		}
		else
		{
			reg = info_bit_padded[_channel_order_sorted.at(nuser).at(j)];
			for (int n = NUM_USER - 1; n >= 0; n--)
			{
				info_bit_binary.at(n).at(j) = reg % 2;
				reg /= 2;
			}
		}
	}
	bool crc_pass = true;
	
	for (int n = 0; n < NUM_USER; n++)
	{
		for (uint16_t i = _info_length; i < _info_length + _crc_size; ++i)
		{
			uint8_t crc_bit = 0;
			for (uint16_t j = 0; j < _info_length; ++j)
			{
				crc_bit = (uint8_t)((crc_bit + _crc_matrix.at(i - _info_length).at(j) * info_bit_binary.at(n).at(j)) % 2);
			}
			
			if (crc_bit != info_bit_binary.at(n).at(i))
			{
				crc_pass = false;
				break;
			}
		}
		if (!JOINT || crc_pass == false)
			break;
	}
	/*for (uint16_t i = _info_length; i < _info_length + _crc_size; ++i)
	{
		uint8_t crc_bit = 0;
		for (uint16_t j = 0; j < _info_length; ++j)
		{
			crc_bit = (uint8_t)((crc_bit + _crc_matrix.at(i - _info_length).at(j) * info_bit_padded[_channel_order_sorted.at(nuser).at(j)]) % 2);
		}
		if (crc_bit != info_bit_padded[_channel_order_sorted.at(nuser).at(i)])
		{
			crc_pass = false;
			break;
		}
	}*/
		
	return crc_pass;
}

std::vector<int> PolarCode::decode_scl_llr(uint16_t list_size, double *llr, int nuser)
{
	for (int i = 0; i < _block_length; i++)
	{
		if (isinf(llr[i]))
		{
			llr[i] = signbit(llr[i]) ? -40 : 40;
		}
		llr[i] > 40 ? llr[i] = 40 : 0;
		llr[i] < -40 ? llr[i] = -40 : 0;
	}
	_list_size = list_size;
	_llr_based_computation = true;
	initializeDataStructures();
	uint16_t  l = assignInitialPath();
	double *llr_0 = getArrayPointer_LLR(0, l);
	for (uint16_t beta = 0; beta < _block_length; ++beta)
	{
		llr_0[beta] = llr[beta];
	}

	return decode_scl(nuser);
}

void PolarCode::decode_jpscl_llr(uint16_t list_size, double** app, std::vector<std::vector<int>>& decodedResult, int **Interleaver)
{
	if (NUM_USER == 1)
	{
		std :: cout << "paramater set in decoded_jpsck_llr";
		//system("pause");
	}
	

	/////////////////////////////////////////////////////////
	std::vector<std::vector<int>> Interleaver_invert(NUM_USER, std::vector<int>(CODE_LEN));
	if (INTERLEAVER)
	{
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			for (int i = 0; i < CODE_LEN; i++)
			{
				Interleaver_invert[nuser][Interleaver[nuser][i]] = i;
			}
		}
	}

	_list_size = list_size;
	_llr_based_computation = true;
	initializeDataStructures();
	uint16_t  l = assignInitialPath();
	double** llr_0 = getArrayPointer_pr_vector(0,l);
	for (uint16_t beta = 0; beta < _block_length; ++beta)
	{
		for (int i = 0; i < NUM_LEVEL; i++)
		{
			llr_0[i][beta] = app[Interleaver_invert[1][beta]][i];
		}
	}

	for (uint16_t phi = 0; phi < _block_length; ++phi)
	{
		if (_llr_based_computation)
		{
			recursivelyCalcLLR(_n, phi);  //�ק�message passing���覡
		}
		else
		{
			recursivelyCalcP(_n, phi);  //�ثe�Τ���
		}
		if (_frozen_bits.at(0).at(phi) == 1)
		{
			continuePaths_FrozenBit(phi);
		}
		else
		{
			continuePaths_UnfrozenBit(phi);
		}
		if ((phi % 2) == 1)
		{
			recursivelyUpdateC(_n, phi);
		}
	}
	l = findMostProbablePath((bool)_crc_size, 1);
	uint8_t* c_0 = _arrayPointer_Info.at(l);
	std::vector<uint8_t> deocded_info_bits(_info_length);
	
	for (uint16_t beta = 0; beta < _info_length; ++beta)
	{
		//std::cout << int(c_0[_channel_order_sorted.at(0).at(beta)]) << " ";
		/*decodedResult[0][beta] = c_0[_channel_order_sorted.at(0).at(beta)] / 4;
		decodedResult[1][beta] = c_0[_channel_order_sorted.at(0).at(beta)] / 2;
		decodedResult[2][beta] = c_0[_channel_order_sorted.at(0).at(beta)] % 2;*/
		int reg = c_0[_channel_order_sorted.at(0).at(beta)];
		for (int nuser = NUM_USER - 1; nuser >= 0; nuser--)
		{
			decodedResult[nuser][beta] = reg % 2;
			reg /= 2;
		}

	}
	//system("pause");
	

	for (uint16_t s = 0; s < _list_size; ++s)
	{
		delete[] _arrayPointer_Info.at(s);
		for (uint16_t lambda = 0; lambda < _n + 1; ++lambda)
		{
			for (int i = 0; i < NUM_LEVEL; i++)
			{
				delete[] _arrayPointer_pr_vector.at(lambda).at(s)[i];
				delete[] _arrayPointer_C_vector.at(lambda).at(s)[i];
			}
			delete[] _arrayPointer_pr_vector.at(lambda).at(s);
			delete[] _arrayPointer_C_vector.at(lambda).at(s);
		}
	}
	
	return ;
}

std::vector<int> PolarCode::decode_scl(int nuser)
{
	for (uint16_t phi = 0; phi < _block_length; ++phi)
	{
		if (_llr_based_computation)
		{
			recursivelyCalcLLR(_n, phi);  //�ק�message passing���覡
		}
		else 
		{
			recursivelyCalcP(_n, phi);  //�ثe�Τ���
		}
		if (NBC && _frozen_bits.at(nuser).at(phi) == 1 && phi != (1023 - 960) )
		{
			continuePaths_FrozenBit(phi);
		}
		else if (!NBC && _frozen_bits.at(nuser).at(phi) == 1)
		{
			continuePaths_FrozenBit(phi);
		}
		else
		{
			continuePaths_UnfrozenBit(phi);
		}
		if ((phi % 2) == 1)
		{
			recursivelyUpdateC(_n, phi);
		}
	}
	
	uint16_t l = findMostProbablePath((bool)_crc_size, nuser);
	uint8_t * c_0 = _arrayPointer_Info.at(l);
	std::vector<int> deocded_info_bits(_info_length);
	for (uint16_t beta = 0; beta < _info_length; ++beta)
	{
		deocded_info_bits[beta] = int(c_0[_channel_order_sorted.at(nuser).at(beta)]);
	}
	for (uint16_t s = 0; s < _list_size; ++s)
	{
		delete[] _arrayPointer_Info.at(s);
		for (uint16_t lambda = 0; lambda < _n + 1; ++lambda)
		{
			if (_llr_based_computation)
			{
				delete[] _arrayPointer_LLR.at(lambda).at(s);
			}
			else
			{
				delete[] _arrayPointer_P.at(lambda).at(s);
			}
			delete[] _arrayPointer_C.at(lambda).at(s);
		}
	}
	return deocded_info_bits;
}

void PolarCode::initializeDataStructures()
{
	while (_inactivePathIndices.size())
	{
		_inactivePathIndices.pop();
	};
	_activePath.resize(_list_size);
	if (!JOINT)
	{
		if (_llr_based_computation)
		{
			_pathMetric_LLR.resize(_list_size);
			_arrayPointer_LLR.resize(_n + 1);
			for (int i = 0; i < _n + 1; ++i)
			{
				_arrayPointer_LLR.at(i).resize(_list_size);
			}
		}
		else
		{
			_arrayPointer_P.resize(_n + 1);
			for (int i = 0; i < _n + 1; ++i)
			{
				_arrayPointer_P.at(i).resize(_list_size);
			}
		}
		_arrayPointer_C.resize(_n + 1);
		for (int i = 0; i < _n + 1; ++i)
		{
			_arrayPointer_C.at(i).resize(_list_size);
		}
	}
	else
	{
		// initial pr_vector and C_vector
		_pathMetric_LLR.resize(_list_size);
		_arrayPointer_pr_vector.resize(_n + 1);
		for (int i = 0; i < _n + 1; ++i)
		{
			_arrayPointer_pr_vector.at(i).resize(_list_size);
			for (int j = 0; j < _list_size; j++)
			{
				_arrayPointer_pr_vector.at(i).at(j) = new double*[NUM_LEVEL];
			}
		}
		_arrayPointer_C_vector.resize(_n + 1);
		for (int i = 0; i < _n + 1; ++i)
		{
			_arrayPointer_C_vector.at(i).resize(_list_size);
			for (int j = 0; j < _list_size; j++)
			{
				_arrayPointer_C_vector.at(i).at(j) = new double* [NUM_LEVEL];
			}
		}
	}
	_arrayPointer_Info.resize(_list_size);
	_pathIndexToArrayIndex.resize(_n + 1);
	for (int i = 0; i < _n + 1; ++i)
	{
		_pathIndexToArrayIndex.at(i).resize(_list_size);
	}
	_inactiveArrayIndices.resize(_n + 1);
	for (int i = 0; i < _n + 1; ++i)
	{
		while (_inactiveArrayIndices.at(i).size())
		{
			_inactiveArrayIndices.at(i).pop();
		}
	}
	_arrayReferenceCount.resize(_n + 1);
	for (int i = 0; i < _n + 1; ++i)
	{
		_arrayReferenceCount.at(i).resize(_list_size);
	}
	for (uint16_t s = 0; s < _list_size; ++s)
	{
		_arrayPointer_Info.at(s) = new uint8_t[_block_length](); //_arrayPointer_Info
		for (uint16_t lambda = 0; lambda < _n + 1; ++lambda)
		{
			if (!JOINT)
			{
				if (_llr_based_computation)
				{
					_arrayPointer_LLR.at(lambda).at(s) = new double[(1 << (_n - lambda))]();
				}
				else
				{
					_arrayPointer_P.at(lambda).at(s) = new double[2 * (1 << (_n - lambda))]();
				}
				_arrayPointer_C.at(lambda).at(s) = new uint8_t[2 * (1 << (_n - lambda))]();
			}
			else
			{
				for (int i = 0; i < NUM_LEVEL; i++)
				{
					_arrayPointer_pr_vector.at(lambda).at(s)[i] = new double[(1 << (_n - lambda))]();
					_arrayPointer_C_vector.at(lambda).at(s)[i] = new double[2 * (1 << (_n - lambda))]();
				}
			}
			
			_arrayReferenceCount.at(lambda).at(s) = 0;
			_inactiveArrayIndices.at(lambda).push(s);
		}
	}
	for (uint16_t l = 0; l < _list_size; ++l)
	{
		_activePath.at(l) = 0;
		_inactivePathIndices.push(l);
		if (_llr_based_computation)
		{
			_pathMetric_LLR.at(l) = 0;
		}
	}
}

uint16_t PolarCode::assignInitialPath()
{
	uint16_t l = _inactivePathIndices.top();
	_inactivePathIndices.pop();
	_activePath.at(l) = 1;
	// Associate arrays with path index
	for (uint16_t lambda = 0; lambda < _n + 1; ++lambda)
	{
		uint16_t  s = _inactiveArrayIndices.at(lambda).top();
		_inactiveArrayIndices.at(lambda).pop();
		_pathIndexToArrayIndex.at(lambda).at(l) = s;
		_arrayReferenceCount.at(lambda).at(s) = 1;
	}
	return l;
}

uint16_t PolarCode::clonePath(uint16_t l)
{
	uint16_t l_p = _inactivePathIndices.top();
	_inactivePathIndices.pop();
	_activePath.at(l_p) = 1;
	if (!JOINT && _llr_based_computation)
	{
		_pathMetric_LLR.at(l_p) = _pathMetric_LLR.at(l);
	}
	for (uint16_t lambda = 0; lambda < _n + 1; ++lambda)
	{
		uint16_t s = _pathIndexToArrayIndex.at(lambda).at(l);
		_pathIndexToArrayIndex.at(lambda).at(l_p) = s;
		_arrayReferenceCount.at(lambda).at(s)++;
	}
	return l_p;
}

void PolarCode::killPath(uint16_t l)
{
	_activePath.at(l) = 0;
	_inactivePathIndices.push(l);
	if (_llr_based_computation)
	{
		_pathMetric_LLR.at(l) = 0;
	}
	for (uint16_t lambda = 0; lambda < _n + 1; ++lambda)
	{
		uint16_t s = _pathIndexToArrayIndex.at(lambda).at(l);
		_arrayReferenceCount.at(lambda).at(s)--;
		if (_arrayReferenceCount.at(lambda).at(s) == 0)
		{
			_inactiveArrayIndices.at(lambda).push(s);
		}
	}
}

double *PolarCode::getArrayPointer_P(uint16_t lambda, uint16_t l)
{
	uint16_t s = _pathIndexToArrayIndex.at(lambda).at(l);
	uint16_t s_p;
	if (_arrayReferenceCount.at(lambda).at(s) == 1)
	{
		s_p = s;
	}
	else
	{
		s_p = _inactiveArrayIndices.at(lambda).top();
		_inactiveArrayIndices.at(lambda).pop();
		//copy
		std::copy(_arrayPointer_P.at(lambda).at(s), _arrayPointer_P.at(lambda).at(s) + (1 << (_n - lambda + 1)), _arrayPointer_P.at(lambda).at(s_p));
		std::copy(_arrayPointer_C.at(lambda).at(s), _arrayPointer_C.at(lambda).at(s) + (1 << (_n - lambda + 1)), _arrayPointer_C.at(lambda).at(s_p));
		_arrayReferenceCount.at(lambda).at(s)--;
		_arrayReferenceCount.at(lambda).at(s_p) = 1;
		_pathIndexToArrayIndex.at(lambda).at(l) = s_p;
	}
	return _arrayPointer_P.at(lambda).at(s_p);
}

double *PolarCode::getArrayPointer_LLR(uint16_t lambda, uint16_t l)
{
	uint16_t s = _pathIndexToArrayIndex.at(lambda).at(l);
	uint16_t s_p;
	if (_arrayReferenceCount.at(lambda).at(s) == 1)
	{
		s_p = s;
	}
	else
	{
		s_p = _inactiveArrayIndices.at(lambda).top();
		_inactiveArrayIndices.at(lambda).pop();
		//copy
		std::copy(_arrayPointer_C.at(lambda).at(s), _arrayPointer_C.at(lambda).at(s) + (1 << (_n - lambda + 1)), _arrayPointer_C.at(lambda).at(s_p));
		std::copy(_arrayPointer_LLR.at(lambda).at(s), _arrayPointer_LLR.at(lambda).at(s) + (1 << (_n - lambda)), _arrayPointer_LLR.at(lambda).at(s_p));

		_arrayReferenceCount.at(lambda).at(s)--;
		_arrayReferenceCount.at(lambda).at(s_p) = 1;
		_pathIndexToArrayIndex.at(lambda).at(l) = s_p;
	}
	return _arrayPointer_LLR.at(lambda).at(s_p);
}

double **PolarCode::getArrayPointer_pr_vector(uint16_t lambda, uint16_t l)
{
	uint16_t s = _pathIndexToArrayIndex.at(lambda).at(l);
	uint16_t s_p;
	if (_arrayReferenceCount.at(lambda).at(s) == 1)
	{
		s_p = s;
	}
	else
	{
		s_p = _inactiveArrayIndices.at(lambda).top();
		_inactiveArrayIndices.at(lambda).pop();
		//----copy
		
		for (int i = 0; i < (1 << (_n - lambda + 1)); i++)
		{
			for (int j = 0; j < NUM_LEVEL; j++)
			{
				_arrayPointer_C_vector.at(lambda).at(s_p)[j][i] = _arrayPointer_C_vector.at(lambda).at(s)[j][i];
			}
		}
		for (int i = 0; i < (1 << (_n - lambda)); i++)
		{
			for (int j = 0; j < NUM_LEVEL; j++)
			{
				_arrayPointer_pr_vector.at(lambda).at(s_p)[j][i] = _arrayPointer_pr_vector.at(lambda).at(s)[j][i];
			}
		}
		_arrayReferenceCount.at(lambda).at(s)--;
		_arrayReferenceCount.at(lambda).at(s_p) = 1;
		_pathIndexToArrayIndex.at(lambda).at(l) = s_p;
	}
	return _arrayPointer_pr_vector.at(lambda).at(s_p);
}

uint8_t *PolarCode::getArrayPointer_C(uint16_t lambda, uint16_t l)
{
	uint16_t s = _pathIndexToArrayIndex.at(lambda).at(l);
	uint16_t s_p;
	if (_arrayReferenceCount.at(lambda).at(s) == 1)
	{
		s_p = s;
	}
	else
	{
		s_p = _inactiveArrayIndices.at(lambda).top();
		_inactiveArrayIndices.at(lambda).pop();
		//copy
		if (_llr_based_computation)
		{
			std::copy(_arrayPointer_LLR.at(lambda).at(s), _arrayPointer_LLR.at(lambda).at(s) + (1 << (_n - lambda)), _arrayPointer_LLR.at(lambda).at(s_p));
		}
		else
		{
			std::copy(_arrayPointer_P.at(lambda).at(s), _arrayPointer_P.at(lambda).at(s) + (1 << (_n - lambda + 1)), _arrayPointer_P.at(lambda).at(s_p));
		}
		std::copy(_arrayPointer_C.at(lambda).at(s), _arrayPointer_C.at(lambda).at(s) + (1 << (_n - lambda + 1)), _arrayPointer_C.at(lambda).at(s_p));
		_arrayReferenceCount.at(lambda).at(s)--;
		_arrayReferenceCount.at(lambda).at(s_p) = 1;
		_pathIndexToArrayIndex.at(lambda).at(l) = s_p;
	}
	return _arrayPointer_C.at(lambda).at(s_p);
}

double **PolarCode::getArrayPointer_C_vector(uint16_t lambda, uint16_t l)
{
	uint16_t s = _pathIndexToArrayIndex.at(lambda).at(l);
	uint16_t s_p;
	if (_arrayReferenceCount.at(lambda).at(s) == 1)
	{
		s_p = s;
	}
	else
	{
		s_p = _inactiveArrayIndices.at(lambda).top();
		_inactiveArrayIndices.at(lambda).pop();
		//----copy
		if (_llr_based_computation)
		{
			for (int i = 0; i < (1 << (_n - lambda)); i++)
			{
				for (int j = 0; j < NUM_LEVEL; j++)
				{
					_arrayPointer_pr_vector.at(lambda).at(s_p)[j][i] = _arrayPointer_pr_vector.at(lambda).at(s)[j][i];
				}
			}
		}
		else
		{
			std::cout << "paramater set error !!! in getArrayPointer_C_vector";
			//system("pause");
		}
		
		for (int i = 0; i < (1 << (_n - lambda + 1)); i++)
		{
			for (int j = 0; j < NUM_LEVEL; j++)
			{
				_arrayPointer_C_vector.at(lambda).at(s_p)[j][i] = _arrayPointer_C_vector.at(lambda).at(s)[j][i];
			}
		}
		_arrayReferenceCount.at(lambda).at(s)--;
		_arrayReferenceCount.at(lambda).at(s_p) = 1;
		_pathIndexToArrayIndex.at(lambda).at(l) = s_p;
	}
	return _arrayPointer_C_vector.at(lambda).at(s_p);
}

void PolarCode::recursivelyCalcP(uint16_t lambda, uint16_t phi)
{
	std::cout << "recursivelyCalcP";
	//system("pause");
	if (lambda == 0)
	{
		return;
	}
	uint16_t psi = phi >> 1;
	if ((phi % 2) == 0)
	{
		recursivelyCalcP(lambda - 1, psi);
	}
	double sigma = 0.0f;
	for (uint16_t l = 0; l < _list_size; ++l)
	{
		if (_activePath.at(l) == 0)
		{
			continue;
		}
		double *p_lambda = getArrayPointer_P(lambda, l);
		double *p_lambda_1 = getArrayPointer_P(lambda - 1, l);
		uint8_t * c_lambda = getArrayPointer_C(lambda, l);
		for (uint16_t beta = 0; beta < (1 << (_n - lambda)); ++beta)
		{
			if ((phi % 2) == 0)
			{
				p_lambda[2 * beta] = 0.5f * (p_lambda_1[2 * (2 * beta)] * p_lambda_1[2 * (2 * beta + 1)]
					+ p_lambda_1[2 * (2 * beta) + 1] * p_lambda_1[2 * (2 * beta + 1) + 1]);
				p_lambda[2 * beta + 1] = 0.5f * (p_lambda_1[2 * (2 * beta) + 1] * p_lambda_1[2 * (2 * beta + 1)]
					+ p_lambda_1[2 * (2 * beta)] * p_lambda_1[2 * (2 * beta + 1) + 1]);
			}
			else
			{
				uint8_t  u_p = c_lambda[2 * beta];
				p_lambda[2 * beta] = 0.5f * p_lambda_1[2 * (2 * beta) + (u_p % 2)] * p_lambda_1[2 * (2 * beta + 1)];
				p_lambda[2 * beta + 1] = 0.5f * p_lambda_1[2 * (2 * beta) + ((u_p + 1) % 2)] * p_lambda_1[2 * (2 * beta + 1) + 1];
			}
			sigma = std::max(sigma, p_lambda[2 * beta]);
			sigma = std::max(sigma, p_lambda[2 * beta + 1]);
		}
	}
	for (uint16_t l = 0; l < _list_size; ++l)
	{
		if (sigma == 0) // Typically happens because of undeflow
		{
			break;
		}
		if (_activePath.at(l) == 0)
		{
			continue;
		}
		double *p_lambda = getArrayPointer_P(lambda, l);
		for (uint16_t beta = 0; beta < (1 << (_n - lambda)); ++beta)
		{
			p_lambda[2 * beta] = p_lambda[2 * beta] / sigma;
			p_lambda[2 * beta + 1] = p_lambda[2 * beta + 1] / sigma;
		}
	}
}

void PolarCode::recursivelyCalcLLR(uint16_t lambda, uint16_t phi)
{
	// lambda is layer 0 ~ BCT_LAYER- 1 , phi is phase 0 ~ CODE_LEN - 1  
	// lambda table polar code (8,4)
	// 0 0 0 0
	// 1 0 0 0
	// 2 1 0 0
	// 3 1 0 0
	// 4 2 1 0
	// 5 2 1 0
	// .......
	if (lambda == 0) // stoping condition
	{
		return;
	}
	uint16_t psi = phi >> 1;
	if ((phi % 2) == 0)
	{
		recursivelyCalcLLR(lambda - 1, psi);
	}
	for (uint16_t l = 0; l < _list_size; ++l)
	{
		if (_activePath.at(l) == 0)
		{
			continue;
		}
		if (!JOINT) //---SCL
		{
			double* llr_lambda = getArrayPointer_LLR(lambda, l);
			double* llr_lambda_1 = getArrayPointer_LLR(lambda - 1, l);
			uint8_t* c_lambda = getArrayPointer_C(lambda, l);
			for (uint16_t beta = 0; beta < (1 << (_n - lambda)); ++beta)
			{
				if ((phi % 2) == 0) // upper branch
				{					// min sum
					if (40 > std::max(std::abs(llr_lambda_1[2 * beta]), std::abs(llr_lambda_1[2 * beta + 1])))
					{
						llr_lambda[beta] = std::log((exp(llr_lambda_1[2 * beta] + llr_lambda_1[2 * beta + 1]) + 1) /
							(exp(llr_lambda_1[2 * beta]) + exp(llr_lambda_1[2 * beta + 1])));
					}
					else
					{
						llr_lambda[beta] = (double)((llr_lambda_1[2 * beta] < 0) ? -1 : (llr_lambda_1[2 * beta] > 0)) *
							((llr_lambda_1[2 * beta + 1] < 0) ? -1 : (llr_lambda_1[2 * beta + 1] > 0)) *
							std::min(std::abs(llr_lambda_1[2 * beta]), std::abs(llr_lambda_1[2 * beta + 1]));
					}
				}
				else // lower branch
				{
					uint8_t u_p = c_lambda[2 * beta];
					llr_lambda[beta] = (1 - 2 * u_p) * llr_lambda_1[2 * beta] + llr_lambda_1[2 * beta + 1];
				}
			}
		}
		else    //----JPSCL
		{
			double** llr_lambda = getArrayPointer_pr_vector(lambda, l);
			double** llr_lambda_1 = getArrayPointer_pr_vector(lambda - 1, l);
			double** c_lambda = getArrayPointer_C_vector(lambda, l);

			for (uint16_t beta = 0; beta < (1 << (_n - lambda)); ++beta)
			{
				if ((phi % 2) == 0) // upper branch
				{
					//---- f function
					if (NUM_LEVEL == 4)
					{
						llr_lambda[0][beta] = llr_lambda_1[0][2 * beta] * llr_lambda_1[0][2 * beta + 1] + llr_lambda_1[1][2 * beta] * llr_lambda_1[1][2 * beta + 1] + llr_lambda_1[2][2 * beta] * llr_lambda_1[2][2 * beta + 1] + llr_lambda_1[3][2 * beta] * llr_lambda_1[3][2 * beta + 1];
						llr_lambda[1][beta] = llr_lambda_1[0][2 * beta] * llr_lambda_1[1][2 * beta + 1] + llr_lambda_1[1][2 * beta] * llr_lambda_1[0][2 * beta + 1] + llr_lambda_1[2][2 * beta] * llr_lambda_1[3][2 * beta + 1] + llr_lambda_1[3][2 * beta] * llr_lambda_1[2][2 * beta + 1];
						llr_lambda[2][beta] = llr_lambda_1[0][2 * beta] * llr_lambda_1[2][2 * beta + 1] + llr_lambda_1[1][2 * beta] * llr_lambda_1[3][2 * beta + 1] + llr_lambda_1[2][2 * beta] * llr_lambda_1[0][2 * beta + 1] + llr_lambda_1[3][2 * beta] * llr_lambda_1[1][2 * beta + 1];
						llr_lambda[3][beta] = llr_lambda_1[0][2 * beta] * llr_lambda_1[3][2 * beta + 1] + llr_lambda_1[1][2 * beta] * llr_lambda_1[2][2 * beta + 1] + llr_lambda_1[2][2 * beta] * llr_lambda_1[1][2 * beta + 1] + llr_lambda_1[3][2 * beta] * llr_lambda_1[0][2 * beta + 1];
					}
					else if (NUM_LEVEL == 8)
					{
						llr_lambda[0][beta] = llr_lambda_1[0][2 * beta] * llr_lambda_1[0][2 * beta + 1] + llr_lambda_1[1][2 * beta] * llr_lambda_1[1][2 * beta + 1] + llr_lambda_1[2][2 * beta] * llr_lambda_1[2][2 * beta + 1] + llr_lambda_1[3][2 * beta] * llr_lambda_1[3][2 * beta + 1] + llr_lambda_1[4][2 * beta] * llr_lambda_1[4][2 * beta + 1] + llr_lambda_1[5][2 * beta] * llr_lambda_1[5][2 * beta + 1] + llr_lambda_1[6][2 * beta] * llr_lambda_1[6][2 * beta + 1] + llr_lambda_1[7][2 * beta] * llr_lambda_1[7][2 * beta + 1];
						llr_lambda[1][beta] = llr_lambda_1[0][2 * beta] * llr_lambda_1[1][2 * beta + 1] + llr_lambda_1[1][2 * beta] * llr_lambda_1[0][2 * beta + 1] + llr_lambda_1[2][2 * beta] * llr_lambda_1[3][2 * beta + 1] + llr_lambda_1[3][2 * beta] * llr_lambda_1[2][2 * beta + 1] + llr_lambda_1[4][2 * beta] * llr_lambda_1[5][2 * beta + 1] + llr_lambda_1[5][2 * beta] * llr_lambda_1[4][2 * beta + 1] + llr_lambda_1[6][2 * beta] * llr_lambda_1[7][2 * beta + 1] + llr_lambda_1[7][2 * beta] * llr_lambda_1[6][2 * beta + 1];
						llr_lambda[2][beta] = llr_lambda_1[0][2 * beta] * llr_lambda_1[2][2 * beta + 1] + llr_lambda_1[1][2 * beta] * llr_lambda_1[3][2 * beta + 1] + llr_lambda_1[2][2 * beta] * llr_lambda_1[0][2 * beta + 1] + llr_lambda_1[3][2 * beta] * llr_lambda_1[1][2 * beta + 1] + llr_lambda_1[4][2 * beta] * llr_lambda_1[6][2 * beta + 1] + llr_lambda_1[5][2 * beta] * llr_lambda_1[7][2 * beta + 1] + llr_lambda_1[6][2 * beta] * llr_lambda_1[4][2 * beta + 1] + llr_lambda_1[7][2 * beta] * llr_lambda_1[5][2 * beta + 1];
						llr_lambda[3][beta] = llr_lambda_1[0][2 * beta] * llr_lambda_1[3][2 * beta + 1] + llr_lambda_1[1][2 * beta] * llr_lambda_1[2][2 * beta + 1] + llr_lambda_1[2][2 * beta] * llr_lambda_1[1][2 * beta + 1] + llr_lambda_1[3][2 * beta] * llr_lambda_1[0][2 * beta + 1] + llr_lambda_1[4][2 * beta] * llr_lambda_1[7][2 * beta + 1] + llr_lambda_1[5][2 * beta] * llr_lambda_1[6][2 * beta + 1] + llr_lambda_1[6][2 * beta] * llr_lambda_1[5][2 * beta + 1] + llr_lambda_1[7][2 * beta] * llr_lambda_1[4][2 * beta + 1];
						llr_lambda[4][beta] = llr_lambda_1[0][2 * beta] * llr_lambda_1[4][2 * beta + 1] + llr_lambda_1[1][2 * beta] * llr_lambda_1[5][2 * beta + 1] + llr_lambda_1[2][2 * beta] * llr_lambda_1[6][2 * beta + 1] + llr_lambda_1[3][2 * beta] * llr_lambda_1[7][2 * beta + 1] + llr_lambda_1[4][2 * beta] * llr_lambda_1[0][2 * beta + 1] + llr_lambda_1[5][2 * beta] * llr_lambda_1[1][2 * beta + 1] + llr_lambda_1[6][2 * beta] * llr_lambda_1[2][2 * beta + 1] + llr_lambda_1[7][2 * beta] * llr_lambda_1[3][2 * beta + 1];
						llr_lambda[5][beta] = llr_lambda_1[0][2 * beta] * llr_lambda_1[5][2 * beta + 1] + llr_lambda_1[1][2 * beta] * llr_lambda_1[4][2 * beta + 1] + llr_lambda_1[2][2 * beta] * llr_lambda_1[7][2 * beta + 1] + llr_lambda_1[3][2 * beta] * llr_lambda_1[6][2 * beta + 1] + llr_lambda_1[4][2 * beta] * llr_lambda_1[1][2 * beta + 1] + llr_lambda_1[5][2 * beta] * llr_lambda_1[0][2 * beta + 1] + llr_lambda_1[6][2 * beta] * llr_lambda_1[3][2 * beta + 1] + llr_lambda_1[7][2 * beta] * llr_lambda_1[2][2 * beta + 1];
						llr_lambda[6][beta] = llr_lambda_1[0][2 * beta] * llr_lambda_1[6][2 * beta + 1] + llr_lambda_1[1][2 * beta] * llr_lambda_1[7][2 * beta + 1] + llr_lambda_1[2][2 * beta] * llr_lambda_1[4][2 * beta + 1] + llr_lambda_1[3][2 * beta] * llr_lambda_1[5][2 * beta + 1] + llr_lambda_1[4][2 * beta] * llr_lambda_1[2][2 * beta + 1] + llr_lambda_1[5][2 * beta] * llr_lambda_1[3][2 * beta + 1] + llr_lambda_1[6][2 * beta] * llr_lambda_1[0][2 * beta + 1] + llr_lambda_1[7][2 * beta] * llr_lambda_1[1][2 * beta + 1];
						llr_lambda[7][beta] = llr_lambda_1[0][2 * beta] * llr_lambda_1[7][2 * beta + 1] + llr_lambda_1[1][2 * beta] * llr_lambda_1[6][2 * beta + 1] + llr_lambda_1[2][2 * beta] * llr_lambda_1[5][2 * beta + 1] + llr_lambda_1[3][2 * beta] * llr_lambda_1[4][2 * beta + 1] + llr_lambda_1[4][2 * beta] * llr_lambda_1[3][2 * beta + 1] + llr_lambda_1[5][2 * beta] * llr_lambda_1[2][2 * beta + 1] + llr_lambda_1[6][2 * beta] * llr_lambda_1[1][2 * beta + 1] + llr_lambda_1[7][2 * beta] * llr_lambda_1[0][2 * beta + 1];
					}
					//----limit bound
					double sum = 0;
					for (int j = 0; j < NUM_LEVEL; j++)
					{
						if (llr_lambda[j][beta] < NUMERIC_LIMIT)
							llr_lambda[j][beta] = NUMERIC_LIMIT;
						sum += llr_lambda[j][beta];
					}
					//----normal
					for (int l = 0; l < NUM_LEVEL; l++)
					{
						llr_lambda[l][beta] /= sum;
					}
				}
				else // lower branch
				{
					std::vector<double> temp(NUM_LEVEL);
					//---- f function
					if (NUM_LEVEL == 4)
					{
						temp[0] = llr_lambda_1[0][2 * beta] * c_lambda[0][2 * beta] + llr_lambda_1[1][2 * beta] * c_lambda[1][2 * beta] + llr_lambda_1[2][2 * beta] * c_lambda[2][2 * beta] + llr_lambda_1[3][2 * beta] * c_lambda[3][2 * beta];
						temp[1] = llr_lambda_1[0][2 * beta] * c_lambda[1][2 * beta] + llr_lambda_1[1][2 * beta] * c_lambda[0][2 * beta] + llr_lambda_1[2][2 * beta] * c_lambda[3][2 * beta] + llr_lambda_1[3][2 * beta] * c_lambda[2][2 * beta];
						temp[2] = llr_lambda_1[0][2 * beta] * c_lambda[2][2 * beta] + llr_lambda_1[1][2 * beta] * c_lambda[3][2 * beta] + llr_lambda_1[2][2 * beta] * c_lambda[0][2 * beta] + llr_lambda_1[3][2 * beta] * c_lambda[1][2 * beta];
						temp[3] = llr_lambda_1[0][2 * beta] * c_lambda[3][2 * beta] + llr_lambda_1[1][2 * beta] * c_lambda[2][2 * beta] + llr_lambda_1[2][2 * beta] * c_lambda[1][2 * beta] + llr_lambda_1[3][2 * beta] * c_lambda[0][2 * beta];
					}
					else if (NUM_LEVEL == 8)
					{
						temp[0] = llr_lambda_1[0][2 * beta] * c_lambda[0][2 * beta] + llr_lambda_1[1][2 * beta] * c_lambda[1][2 * beta] + llr_lambda_1[2][2 * beta] * c_lambda[2][2 * beta] + llr_lambda_1[3][2 * beta] * c_lambda[3][2 * beta] + llr_lambda_1[4][2 * beta] * c_lambda[4][2 * beta] + llr_lambda_1[5][2 * beta] * c_lambda[5][2 * beta] + llr_lambda_1[6][2 * beta] * c_lambda[6][2 * beta] + llr_lambda_1[7][2 * beta] * c_lambda[7][2 * beta];
						temp[1] = llr_lambda_1[0][2 * beta] * c_lambda[1][2 * beta] + llr_lambda_1[1][2 * beta] * c_lambda[0][2 * beta] + llr_lambda_1[2][2 * beta] * c_lambda[3][2 * beta] + llr_lambda_1[3][2 * beta] * c_lambda[2][2 * beta] + llr_lambda_1[4][2 * beta] * c_lambda[5][2 * beta] + llr_lambda_1[5][2 * beta] * c_lambda[4][2 * beta] + llr_lambda_1[6][2 * beta] * c_lambda[7][2 * beta] + llr_lambda_1[7][2 * beta] * c_lambda[6][2 * beta];
						temp[2] = llr_lambda_1[0][2 * beta] * c_lambda[2][2 * beta] + llr_lambda_1[1][2 * beta] * c_lambda[3][2 * beta] + llr_lambda_1[2][2 * beta] * c_lambda[0][2 * beta] + llr_lambda_1[3][2 * beta] * c_lambda[1][2 * beta] + llr_lambda_1[4][2 * beta] * c_lambda[6][2 * beta] + llr_lambda_1[5][2 * beta] * c_lambda[7][2 * beta] + llr_lambda_1[6][2 * beta] * c_lambda[4][2 * beta] + llr_lambda_1[7][2 * beta] * c_lambda[5][2 * beta];
						temp[3] = llr_lambda_1[0][2 * beta] * c_lambda[3][2 * beta] + llr_lambda_1[1][2 * beta] * c_lambda[2][2 * beta] + llr_lambda_1[2][2 * beta] * c_lambda[1][2 * beta] + llr_lambda_1[3][2 * beta] * c_lambda[0][2 * beta] + llr_lambda_1[4][2 * beta] * c_lambda[7][2 * beta] + llr_lambda_1[5][2 * beta] * c_lambda[6][2 * beta] + llr_lambda_1[6][2 * beta] * c_lambda[5][2 * beta] + llr_lambda_1[7][2 * beta] * c_lambda[4][2 * beta];
						temp[4] = llr_lambda_1[0][2 * beta] * c_lambda[4][2 * beta] + llr_lambda_1[1][2 * beta] * c_lambda[5][2 * beta] + llr_lambda_1[2][2 * beta] * c_lambda[6][2 * beta] + llr_lambda_1[3][2 * beta] * c_lambda[7][2 * beta] + llr_lambda_1[4][2 * beta] * c_lambda[0][2 * beta] + llr_lambda_1[5][2 * beta] * c_lambda[1][2 * beta] + llr_lambda_1[6][2 * beta] * c_lambda[2][2 * beta] + llr_lambda_1[7][2 * beta] * c_lambda[3][2 * beta];
						temp[5] = llr_lambda_1[0][2 * beta] * c_lambda[5][2 * beta] + llr_lambda_1[1][2 * beta] * c_lambda[4][2 * beta] + llr_lambda_1[2][2 * beta] * c_lambda[7][2 * beta] + llr_lambda_1[3][2 * beta] * c_lambda[6][2 * beta] + llr_lambda_1[4][2 * beta] * c_lambda[1][2 * beta] + llr_lambda_1[5][2 * beta] * c_lambda[0][2 * beta] + llr_lambda_1[6][2 * beta] * c_lambda[3][2 * beta] + llr_lambda_1[7][2 * beta] * c_lambda[2][2 * beta];
						temp[6] = llr_lambda_1[0][2 * beta] * c_lambda[6][2 * beta] + llr_lambda_1[1][2 * beta] * c_lambda[7][2 * beta] + llr_lambda_1[2][2 * beta] * c_lambda[4][2 * beta] + llr_lambda_1[3][2 * beta] * c_lambda[5][2 * beta] + llr_lambda_1[4][2 * beta] * c_lambda[2][2 * beta] + llr_lambda_1[5][2 * beta] * c_lambda[3][2 * beta] + llr_lambda_1[6][2 * beta] * c_lambda[0][2 * beta] + llr_lambda_1[7][2 * beta] * c_lambda[1][2 * beta];
						temp[7] = llr_lambda_1[0][2 * beta] * c_lambda[7][2 * beta] + llr_lambda_1[1][2 * beta] * c_lambda[6][2 * beta] + llr_lambda_1[2][2 * beta] * c_lambda[5][2 * beta] + llr_lambda_1[3][2 * beta] * c_lambda[4][2 * beta] + llr_lambda_1[4][2 * beta] * c_lambda[3][2 * beta] + llr_lambda_1[5][2 * beta] * c_lambda[2][2 * beta] + llr_lambda_1[6][2 * beta] * c_lambda[1][2 * beta] + llr_lambda_1[7][2 * beta] * c_lambda[0][2 * beta];
					}
					//----limit bound
					for (int j = 0; j < NUM_LEVEL; j++)
					{
						if (temp[j] < NUMERIC_LIMIT)
							temp[j] = NUMERIC_LIMIT;
					}

					//----g function
					double sum = 0;
					for (int l = 0; l < NUM_LEVEL; l++)
					{
						llr_lambda[l][beta] = temp[l] * llr_lambda_1[l][2 * beta + 1];
						if (llr_lambda[l][beta] < NUMERIC_LIMIT)
							llr_lambda[l][beta] = NUMERIC_LIMIT;

						sum += llr_lambda[l][beta];
					}
					//----normal
					for (int l = 0; l < NUM_LEVEL; l++)
					{
						llr_lambda[l][beta] /= sum;
					}
				}
			}

		}
		
	}
}

void PolarCode::recursivelyUpdateC(uint16_t lambda, uint16_t phi)
{
	uint16_t psi = phi >> 1;
	for (uint16_t l = 0; l < _list_size; ++l)
	{
		if (_activePath.at(l) == 0)
		{
			continue;
		}
		if (!JOINT)
		{
			uint8_t* c_lambda = getArrayPointer_C(lambda, l);
			uint8_t* c_lambda_1 = getArrayPointer_C(lambda - 1, l);
			for (uint16_t beta = 0; beta < (1 << (_n - lambda)); ++beta)
			{
				c_lambda_1[2 * (2 * beta) + (psi % 2)] = (uint8_t)((c_lambda[2 * beta] + c_lambda[2 * beta + 1]) % 2);
				c_lambda_1[2 * (2 * beta + 1) + (psi % 2)] = c_lambda[2 * beta + 1];
			}
		}
		else
		{
			double** c_lambda = getArrayPointer_C_vector(lambda, l);
			double** c_lambda_1 = getArrayPointer_C_vector(lambda - 1, l);
			
			for (uint16_t beta = 0; beta < (1 << (_n - lambda)); ++beta)
			{
				// feedback message vector upper branch (left to right)
				//---- f function
				if (NUM_LEVEL == 4)
				{
					c_lambda_1[0][2 * (2 * beta) + (psi % 2)] = c_lambda[0][2 * beta] * c_lambda[0][2 * beta + 1] + c_lambda[1][2 * beta] * c_lambda[1][2 * beta + 1] + c_lambda[2][2 * beta] * c_lambda[2][2 * beta + 1] + c_lambda[3][2 * beta] * c_lambda[3][2 * beta + 1];
					c_lambda_1[1][2 * (2 * beta) + (psi % 2)] = c_lambda[0][2 * beta] * c_lambda[1][2 * beta + 1] + c_lambda[1][2 * beta] * c_lambda[0][2 * beta + 1] + c_lambda[2][2 * beta] * c_lambda[3][2 * beta + 1] + c_lambda[3][2 * beta] * c_lambda[2][2 * beta + 1];
					c_lambda_1[2][2 * (2 * beta) + (psi % 2)] = c_lambda[0][2 * beta] * c_lambda[2][2 * beta + 1] + c_lambda[1][2 * beta] * c_lambda[3][2 * beta + 1] + c_lambda[2][2 * beta] * c_lambda[0][2 * beta + 1] + c_lambda[3][2 * beta] * c_lambda[1][2 * beta + 1];
					c_lambda_1[3][2 * (2 * beta) + (psi % 2)] = c_lambda[0][2 * beta] * c_lambda[3][2 * beta + 1] + c_lambda[1][2 * beta] * c_lambda[2][2 * beta + 1] + c_lambda[2][2 * beta] * c_lambda[1][2 * beta + 1] + c_lambda[3][2 * beta] * c_lambda[0][2 * beta + 1];
					
				}
				else if (NUM_LEVEL == 8)
				{
					c_lambda_1[0][2 * (2 * beta) + (psi % 2)] = c_lambda[0][2 * beta] * c_lambda[0][2 * beta + 1] + c_lambda[1][2 * beta] * c_lambda[1][2 * beta + 1] + c_lambda[2][2 * beta] * c_lambda[2][2 * beta + 1] + c_lambda[3][2 * beta] * c_lambda[3][2 * beta + 1] + c_lambda[4][2 * beta] * c_lambda[4][2 * beta + 1] + c_lambda[5][2 * beta] * c_lambda[5][2 * beta + 1] + c_lambda[6][2 * beta] * c_lambda[6][2 * beta + 1] + c_lambda[7][2 * beta] * c_lambda[7][2 * beta + 1];
					c_lambda_1[1][2 * (2 * beta) + (psi % 2)] = c_lambda[0][2 * beta] * c_lambda[1][2 * beta + 1] + c_lambda[1][2 * beta] * c_lambda[0][2 * beta + 1] + c_lambda[2][2 * beta] * c_lambda[3][2 * beta + 1] + c_lambda[3][2 * beta] * c_lambda[2][2 * beta + 1] + c_lambda[4][2 * beta] * c_lambda[5][2 * beta + 1] + c_lambda[5][2 * beta] * c_lambda[4][2 * beta + 1] + c_lambda[6][2 * beta] * c_lambda[7][2 * beta + 1] + c_lambda[7][2 * beta] * c_lambda[6][2 * beta + 1];
					c_lambda_1[2][2 * (2 * beta) + (psi % 2)] = c_lambda[0][2 * beta] * c_lambda[2][2 * beta + 1] + c_lambda[1][2 * beta] * c_lambda[3][2 * beta + 1] + c_lambda[2][2 * beta] * c_lambda[0][2 * beta + 1] + c_lambda[3][2 * beta] * c_lambda[1][2 * beta + 1] + c_lambda[4][2 * beta] * c_lambda[6][2 * beta + 1] + c_lambda[5][2 * beta] * c_lambda[7][2 * beta + 1] + c_lambda[6][2 * beta] * c_lambda[4][2 * beta + 1] + c_lambda[7][2 * beta] * c_lambda[5][2 * beta + 1];
					c_lambda_1[3][2 * (2 * beta) + (psi % 2)] = c_lambda[0][2 * beta] * c_lambda[3][2 * beta + 1] + c_lambda[1][2 * beta] * c_lambda[2][2 * beta + 1] + c_lambda[2][2 * beta] * c_lambda[1][2 * beta + 1] + c_lambda[3][2 * beta] * c_lambda[0][2 * beta + 1] + c_lambda[4][2 * beta] * c_lambda[7][2 * beta + 1] + c_lambda[5][2 * beta] * c_lambda[6][2 * beta + 1] + c_lambda[6][2 * beta] * c_lambda[5][2 * beta + 1] + c_lambda[7][2 * beta] * c_lambda[4][2 * beta + 1];
					c_lambda_1[4][2 * (2 * beta) + (psi % 2)] = c_lambda[0][2 * beta] * c_lambda[4][2 * beta + 1] + c_lambda[1][2 * beta] * c_lambda[5][2 * beta + 1] + c_lambda[2][2 * beta] * c_lambda[6][2 * beta + 1] + c_lambda[3][2 * beta] * c_lambda[7][2 * beta + 1] + c_lambda[4][2 * beta] * c_lambda[0][2 * beta + 1] + c_lambda[5][2 * beta] * c_lambda[1][2 * beta + 1] + c_lambda[6][2 * beta] * c_lambda[2][2 * beta + 1] + c_lambda[7][2 * beta] * c_lambda[3][2 * beta + 1];
					c_lambda_1[5][2 * (2 * beta) + (psi % 2)] = c_lambda[0][2 * beta] * c_lambda[5][2 * beta + 1] + c_lambda[1][2 * beta] * c_lambda[4][2 * beta + 1] + c_lambda[2][2 * beta] * c_lambda[7][2 * beta + 1] + c_lambda[3][2 * beta] * c_lambda[6][2 * beta + 1] + c_lambda[4][2 * beta] * c_lambda[1][2 * beta + 1] + c_lambda[5][2 * beta] * c_lambda[0][2 * beta + 1] + c_lambda[6][2 * beta] * c_lambda[3][2 * beta + 1] + c_lambda[7][2 * beta] * c_lambda[2][2 * beta + 1];
					c_lambda_1[6][2 * (2 * beta) + (psi % 2)] = c_lambda[0][2 * beta] * c_lambda[6][2 * beta + 1] + c_lambda[1][2 * beta] * c_lambda[7][2 * beta + 1] + c_lambda[2][2 * beta] * c_lambda[4][2 * beta + 1] + c_lambda[3][2 * beta] * c_lambda[5][2 * beta + 1] + c_lambda[4][2 * beta] * c_lambda[2][2 * beta + 1] + c_lambda[5][2 * beta] * c_lambda[3][2 * beta + 1] + c_lambda[6][2 * beta] * c_lambda[0][2 * beta + 1] + c_lambda[7][2 * beta] * c_lambda[1][2 * beta + 1];
					c_lambda_1[7][2 * (2 * beta) + (psi % 2)] = c_lambda[0][2 * beta] * c_lambda[7][2 * beta + 1] + c_lambda[1][2 * beta] * c_lambda[6][2 * beta + 1] + c_lambda[2][2 * beta] * c_lambda[5][2 * beta + 1] + c_lambda[3][2 * beta] * c_lambda[4][2 * beta + 1] + c_lambda[4][2 * beta] * c_lambda[3][2 * beta + 1] + c_lambda[5][2 * beta] * c_lambda[2][2 * beta + 1] + c_lambda[6][2 * beta] * c_lambda[1][2 * beta + 1] + c_lambda[7][2 * beta] * c_lambda[0][2 * beta + 1];
				}
				/*for (int i = 0; i < 4; i++)
				{
					std::cout << c_lambda[i][2 * (2 * beta) + (psi % 2)] << " ";
				}
				system("pause");*/
				//----limit bound
				double sum = 0;
				for (int j = 0; j < NUM_LEVEL; j++)
				{
					if (c_lambda_1[j][2 * (2 * beta) + (psi % 2)] < NUMERIC_LIMIT)
						c_lambda_1[j][2 * (2 * beta) + (psi % 2)] = NUMERIC_LIMIT;
					sum += c_lambda_1[j][2 * (2 * beta) + (psi % 2)];
				}
				//---- normalization
				for (int j = 0; j < NUM_LEVEL; j++)
				{
					c_lambda_1[j][2 * (2 * beta) + (psi % 2)] /= sum;
				}

				//feedback message vector lower branch (left to right)
				sum = 0;
				for (int j = 0; j < NUM_LEVEL; j++)
				{
					c_lambda_1[j][2 * (2 * beta + 1) + (psi % 2)] = c_lambda[j][2 * beta + 1];
					if (c_lambda_1[j][2 * (2 * beta + 1) + (psi % 2)] < NUMERIC_LIMIT)
						c_lambda_1[j][2 * (2 * beta + 1) + (psi % 2)] = NUMERIC_LIMIT;
					sum += c_lambda_1[j][2 * (2 * beta + 1) + (psi % 2)];
				}

				//---- normalization
				for (int j = 0; j < NUM_LEVEL; j++)
				{
					c_lambda_1[j][2 * (2 * beta + 1) + (psi % 2)] /= sum;
				}
			}
			//system("pause");
		}
	}
	if ((psi % 2) == 1)
	{
		recursivelyUpdateC((uint16_t)(lambda - 1), psi);
	}
}

void PolarCode::continuePaths_FrozenBit(uint16_t phi)
{
	for (uint16_t l = 0; l < _list_size; ++l)
	{
		if (_activePath.at(l) == 0)
		{
			continue;
		}

		if (!JOINT)//---- SCL
		{
			uint8_t* c_m = getArrayPointer_C(_n, l);
			//---- frozen value assumed to be zero
			c_m[(phi % 2)] = 0; 
			
			if (_llr_based_computation)
			{
				double* llr_p = getArrayPointer_LLR(_n, l);
				_pathMetric_LLR.at(l) += log(1 + exp(-llr_p[0]));
			}
		}
		else//---- JPSCL
		{
			double** c_m = getArrayPointer_C_vector(_n, l);

			//---- frozen value assumed to be zero
			//---- �o�̰��]user�� forzen bit�@��
			c_m[0][(phi % 2)] = 1;
			for (int i = 1; i < NUM_LEVEL; i++)
			{
				c_m[i][(phi % 2)] = NUMERIC_LIMIT;
			}
			

			if (_llr_based_computation)
			{
				double** pr_vector = getArrayPointer_pr_vector(_n, l);

				//---- calculate D_N of  each user
				double D_N[NUM_USER] = { 0 };

				for (int nuser = 0; nuser < NUM_USER; nuser++)
				{
					int range = NUM_LEVEL / pow(2, nuser + 1);
					double app_sum[2] = { 0 };
					for (int j = 0; j < NUM_LEVEL; j++)
					{
						if (j % (2 * range) < range)
							app_sum[0] += pr_vector[j][0];
						else
							app_sum[1] += pr_vector[j][0];
					}
					if (app_sum[0] <= NUMERIC_LIMIT) D_N[nuser] = -LLR_LIMIT;
					else if (app_sum[1] <= NUMERIC_LIMIT) D_N[nuser] = LLR_LIMIT;
					else D_N[nuser] = log(app_sum[0] / app_sum[1]);
					
					//---- path metric values for l-th path
					_pathMetric_LLR.at(l) += log(1 + exp(-D_N[nuser]));
				}
			}
		}
		_arrayPointer_Info.at(l)[phi] = 0;
	}
}

void PolarCode::continuePaths_UnfrozenBit(uint16_t phi) 
{
	int fork_size = JOINT ? NUM_LEVEL : 2; //---- fork size of List decoding
	std::vector<double> probForks((unsigned long)(fork_size * _list_size));
	std::vector<uint8_t> contForks((unsigned long)(fork_size * _list_size));
	std::vector<double> probabilities;

	uint16_t i = 0; //---- current list size

	for (unsigned l = 0; l < _list_size; ++l)
	{
		
		if (_activePath.at(l) == 0)
		{
			for(int k = 0 ; k < fork_size ; k++)
				probForks.at(fork_size * l + k) = NAN;
		}
		else
		{
			if (_llr_based_computation)
			{
				if (!JOINT) //---- SCL
				{
					double* llr_p = getArrayPointer_LLR(_n, l);
					if (phi != (1023 - 960))
					{
						probForks.at(2 * l) = -(_pathMetric_LLR.at(l) + log(1 + exp(-llr_p[0])));
						probForks.at(2 * l + 1) = -(_pathMetric_LLR.at(l) + log(1 + exp(llr_p[0])));
					}
					else
					{
						probForks.at(2 * l) = -(_pathMetric_LLR.at(l));
						probForks.at(2 * l + 1) = -(_pathMetric_LLR.at(l));
					}
				}
				else //---- JPSCL
				{
					//---- D_N calculation

					double** pr_vector = getArrayPointer_pr_vector(_n, l);
					double D_N[NUM_USER] = { 0 };
					for (int nuser = 0; nuser < NUM_USER; nuser++)
					{
						int range = NUM_LEVEL / pow(2, nuser + 1);
						double app_sum[2] = { 0 };

						for (int j = 0; j < NUM_LEVEL; j++)
						{
							if (j % (2 * range) < range)
								app_sum[0] += pr_vector[j][0];
							else
								app_sum[1] += pr_vector[j][0];
						}
						if (app_sum[0] <= NUMERIC_LIMIT) D_N[nuser] = -LLR_LIMIT;
						else if (app_sum[1] <= NUMERIC_LIMIT) D_N[nuser] = LLR_LIMIT;
						else D_N[nuser] = log(app_sum[0] / app_sum[1]);
					}
					
					/*probForks.at(4 * l) = -(_pathMetric_LLR.at(l) + log(1 + exp(-D_N[0])) + log(1 + exp(-D_N[1])));
					probForks.at(4 * l + 1) = -(_pathMetric_LLR.at(l) + log(1 + exp(-D_N[0])) + log(1 + exp(D_N[1])));
					probForks.at(4 * l + 2) = -(_pathMetric_LLR.at(l) + log(1 + exp(D_N[0])) + log(1 + exp(-D_N[1])));
					probForks.at(4 * l + 3) = -(_pathMetric_LLR.at(l) + log(1 + exp(D_N[0])) + log(1 + exp(D_N[1])));*/
					
					//---- path forks 
					std::vector<int> binary(NUM_USER);
					for (int i = 0; i < fork_size; i++)
					{
						//---- binary transform
						int reg = i;
						for (int nuser = NUM_USER - 1; nuser >= 0; nuser--)
						{
							binary[nuser] = reg % 2;
							reg /= 2;
						}
						
						//---- add D_N value to each fork, there are "forksize" forks
						for (int nuser = 0; nuser < NUM_USER; nuser++)
							probForks.at(fork_size * l + i) += log(1 + exp(-D_N[nuser] * pow(-1, binary[nuser])));

						probForks.at(fork_size * l + i) += _pathMetric_LLR.at(l);
						probForks.at(fork_size * l + i) *= -1;
					}
				}
			}
			else
			{
				// �o�̥Τ���
				double *p_m = getArrayPointer_P(_n, l);
				probForks.at(2 * l) = p_m[0];
				probForks.at(2 * l + 1) = p_m[1];
			}

			//---- "probabilities" are used for selecting good forks
			for (int k = 0; k < fork_size; k++)
				probabilities.push_back(probForks.at(fork_size * l + k));
			i++;
		}
	}
	

	uint16_t rho = _list_size;
	if ((fork_size * i) < _list_size)
	{
		rho = (uint16_t)fork_size * i;
	}
	for (uint8_t l = 0; l < fork_size * _list_size; ++l)
	{
		contForks.at(l) = 0;
	}
	//---- set a threshold
	std::sort(probabilities.begin(), probabilities.end(), std::greater<double>());
	double threshold = probabilities.at((unsigned long)(rho - 1));

	uint16_t num_paths_continued = 0;
	
	for (uint8_t l = 0; l < fork_size * _list_size; ++l)
	{
		if (probForks.at(l) > threshold)
		{
			contForks.at(l) = 1;
			num_paths_continued++;
		}
		if (num_paths_continued == rho)
		{
			break;
		}
	}
	
	if (num_paths_continued < rho)
	{
		for (uint8_t l = 0; l < fork_size * _list_size; ++l)
		{
			if (probForks.at(l) == threshold)
			{
				contForks.at(l) = 1;
				num_paths_continued++;
			}
			if (num_paths_continued == rho)
			{
				break;
			}
		}
	}

	for (unsigned l = 0; l < _list_size; ++l)
	{
		if (_activePath.at(l) == 0)
		{
			continue;
		}
		if (!JOINT && contForks.at(2 * l) == 0 && contForks.at(2 * l + 1) == 0) //�Q�^�O��path
		{
			killPath(l);
		}
		if (JOINT)
		{
			bool tr = 0;
			for (int i = 0; i < fork_size; i++)
			{
				if (contForks.at(fork_size * l + i) != 0)
				{
					tr = 1;
					break;
				}
			}
			
			if (tr == 1)
				continue;

			killPath(l);
		}
	}

	for (unsigned l = 0; l < _list_size; ++l)
	{
		if (!JOINT) //---- SCL
		{
			if (contForks.at(2 * l) == 0 && contForks.at(2 * l + 1) == 0) //���child���Q�^�O��path
			{
				continue;
			}
			uint8_t* c_m = getArrayPointer_C(_n, l); // �ΨӬ����̫�@�h�� 0 1
			if (contForks.at(2 * l) == 1 && contForks.at(2 * l + 1) == 1) //���child���S�Q�^�O
			{
				c_m[(phi % 2)] = 0;
				uint16_t l_p = clonePath(l);
				c_m = getArrayPointer_C(_n, l_p);
				c_m[(phi % 2)] = 1;
				std::copy(_arrayPointer_Info.at(l), _arrayPointer_Info.at(l) + phi, _arrayPointer_Info.at(l_p)); // _arrayPointer_Info.at[...]���ΨӰO����l-pa
				_arrayPointer_Info.at(l)[phi] = 0;
				_arrayPointer_Info.at(l_p)[phi] = 1;
				if (_llr_based_computation)
				{
					double* llr_p = getArrayPointer_LLR(_n, l);
					_pathMetric_LLR.at(l) += log(1 + exp(-llr_p[0]));
					llr_p = getArrayPointer_LLR(_n, l_p);
					_pathMetric_LLR.at(l_p) += log(1 + exp(llr_p[0]));
				}
			}
			else// �䤤�@��child�Q�^�O
			{
				if (contForks.at(2 * l) == 1)// �Ĥ@��child�Q�^�O�A�]�N�Obit=1�Q�^�O
				{
					c_m[(phi % 2)] = 0;
					_arrayPointer_Info.at(l)[phi] = 0;
					if (_llr_based_computation)
					{
						double* llr_p = getArrayPointer_LLR(_n, l);
						_pathMetric_LLR.at(l) += log(1 + exp(-llr_p[0]));
					}
				}
				else// �ĤG��child�Q�^�O�A�]�N�Obit=0�Q�^�O
				{
					c_m[(phi % 2)] = 1;
					_arrayPointer_Info.at(l)[phi] = 1;
					if (_llr_based_computation)
					{
						double* llr_p = getArrayPointer_LLR(_n, l);
						_pathMetric_LLR.at(l) += log(1 + exp(llr_p[0]));
					}
				}
			}
		}
		else  //---- JPSCL
		{

			//---- if all forks are fail, continue to next list
			bool tr = 0;
			for (int i = 0; i < fork_size; i++)
			{
				if (contForks.at(fork_size * l + i) != 0)
				{
					tr = 1;
					break;
				}
			}
			if (tr == 0)
				continue;

			tr = 0;
			uint16_t l_p = l;
			double temp = _pathMetric_LLR.at(l);
			double** c_m = getArrayPointer_C_vector(_n, l_p); // �ΨӬ����̫�@�h�� �|�[�T�����v
			for (int k = 0; k < fork_size; k++)
			{
				if (contForks.at(fork_size * l + k) == 1) //�P�_4��child�O�_�ŦX���
				{
					if (tr == 1) //�Y�o��path�e�����H�ιL �y�s���@��
					{	
						l_p = clonePath(l); //�y�s��
						_pathMetric_LLR.at(l_p) = temp;
						c_m = getArrayPointer_C_vector(_n, l_p);
						std::copy(_arrayPointer_Info.at(l), _arrayPointer_Info.at(l) + phi, _arrayPointer_Info.at(l_p));//�ƻs�P�@����ƫe�����T��
					}

					for (int j = 0; j < NUM_LEVEL; j++)
					{
						c_m[j][(phi % 2)] = NUMERIC_LIMIT;
					}
					c_m[k][(phi % 2)] = 1;
					_arrayPointer_Info.at(l_p)[phi] = k;
					

					if (_llr_based_computation)
					{
						double** pr_vector = getArrayPointer_pr_vector(_n, l_p);
						//---- D_N calculation
						double D_N[NUM_USER] = { 0 };
						for (int nuser = 0; nuser < NUM_USER; nuser++)
						{
							int range = NUM_LEVEL / pow(2, nuser + 1);
							double app_sum[2] = { 0 };

							for (int j = 0; j < NUM_LEVEL; j++)
							{
								if (j % (2 * range) < range)
									app_sum[0] += pr_vector[j][0];
								else
									app_sum[1] += pr_vector[j][0];
							}
							if (app_sum[0] <= NUMERIC_LIMIT) D_N[nuser] = -LLR_LIMIT;
							else if (app_sum[1] <= NUMERIC_LIMIT) D_N[nuser] = LLR_LIMIT;
							else D_N[nuser] = log(app_sum[0] / app_sum[1]);
						}

						/*if (k == 0)
						{
							_pathMetric_LLR.at(l_p) += log(1 + exp(-D_N[0])) + log(1 + exp(-D_N[1]));
						}
						else if (k == 1)
						{
							_pathMetric_LLR.at(l_p) += log(1 + exp(-D_N[0])) + log(1 + exp(D_N[1]));
						}
						else if (k == 2)
						{
							_pathMetric_LLR.at(l_p) += log(1 + exp(D_N[0])) + log(1 + exp(-D_N[1]));
						}
						else if (k == 3)
						{
							_pathMetric_LLR.at(l_p) += log(1 + exp(D_N[0])) + log(1 + exp(D_N[1]));
						}*/

						
					
						//---- binary transform
						std::vector<int> binary(NUM_USER);
						int reg = k;
						for (int nuser = NUM_USER - 1; nuser >= 0; nuser--)
						{
							binary[nuser] = reg % 2;
							reg /= 2;
						}

						//---- add D_N value to each fork, there are "forksize" forks
						for (int nuser = 0; nuser < NUM_USER; nuser++)
							_pathMetric_LLR.at(l_p) += log(1 + exp(-D_N[nuser] * pow(-1, binary[nuser])));
							
					}
					tr = 1;
				}
				
			}	
		}
	}
}

uint16_t PolarCode::findMostProbablePath(bool check_crc, int nuser)
{
	uint16_t l_p = 0;
	double p_p1 = 0;
	double p_llr = std::numeric_limits<double>::max();
	bool path_with_crc_pass = false;
	for (uint16_t l = 0; l < _list_size; ++l)
	{
		if (_activePath.at(l) == 0)
		{
			continue;
		}
		if ((check_crc) && (!crc_check(_arrayPointer_Info.at(l), nuser)))
		{
		//	std::cout << "CRC??";
		//	system("pause");
			continue;
		}
		path_with_crc_pass = true;
		if (_llr_based_computation)
		{
			if (_pathMetric_LLR.at(l) < p_llr)
			{
				p_llr = _pathMetric_LLR.at(l);
				l_p = l;
			}
		}
		else// �ثe�Τ���
		{
			uint8_t *c_m = getArrayPointer_C(_n, l);
			double *p_m = getArrayPointer_P(_n, l);
			if (p_p1 < p_m[c_m[1]])
			{
				l_p = l;
				p_p1 = p_m[c_m[1]];
			}
		}
	}
	if (path_with_crc_pass)
	{
		return l_p;
	}
	else
	{
		return findMostProbablePath(false, nuser);
	}
}

void PolarCode::create_bit_rev_order()
{
	for (uint16_t i = 0; i < _block_length; ++i)
	{
		uint16_t to_be_reversed = i;
		_bit_rev_order.at(i) = (uint16_t)((to_be_reversed & 1) << (_n - 1));
		for (uint8_t j = (uint8_t)(_n - 1); j; --j)
		{
			to_be_reversed >>= 1;
			_bit_rev_order.at(i) += (to_be_reversed & 1) << (j - 1);
		}
	}
}

double PolarCode::JFunction(double sigma)
{
	return pow((1 - pow(2, -0.3073*pow(sigma, 2 * 0.8935))), 1.1064);
}

double PolarCode::InverseJFunction(double I)
{
	return pow(-(1. / 0.3073)*log2(1 - pow(I, 1. / 1.1064)), 0.5 / 0.8935);
}

void PolarCode::decode_BP_Joint(double** app, std::vector<std::vector<int>>& decodedResult, int ** Interleaver)
{
	std::vector<std::vector<int>> Interleaver_invert(NUM_USER, std::vector<int>(CODE_LEN));
	if (INTERLEAVER)
	{
		for (int nuser = 0; nuser < NUM_USER; nuser++)
		{
			for (int i = 0; i < CODE_LEN; i++)
			{
				Interleaver_invert[nuser][Interleaver[nuser][i]] = i;
			}
		}
	}

	std::vector<std::vector<double>> appLlr(NUM_USER, std::vector<double>(DATA_LEN));
	std::vector<std::vector<std::vector<double>>> R(BCT_LAYER + 1, std::vector<std::vector<double>>(NUM_LEVEL, std::vector<double>(CODE_LEN)));
	std::vector<std::vector<std::vector<double>>> L(BCT_LAYER + 1, std::vector<std::vector<double>>(NUM_LEVEL, std::vector<double>(CODE_LEN)));

	/*for (int i = 0; i < CODE_LEN; i++)
	{
		for (int j = 0; j < NUM_LEVEL; j++)
		{
			std::cout << app[i][j] << " ";
		}
		std::cout << std::endl;
	}*/

	//---Initailization
	for (int i = 0; i < CODE_LEN; i++)
	{
		for (int j = 0; j < BCT_LAYER + 1; j++)
		{
			for (int k = 0; k < NUM_LEVEL; k++)
			{
				L[j][k][i] = 1 / double(NUM_LEVEL);
				if (_frozen_bits[0][i] && j == BCT_LAYER)
				{
					R[BCT_LAYER][0][i] = 1;
					for (int l = 1; l < NUM_LEVEL; l++)
						R[BCT_LAYER][l][i] = NUMERIC_LIMIT;
					break;
				}
				else
				{
					R[j][k][i] = 1 / double(NUM_LEVEL);
				}
			}
		}

		for (int j = 0; j < NUM_LEVEL; j++)
			L[0][j][i] = app[Interleaver_invert[1][i]][j];

	}



	/*	for (int i = 0; i < CODE_LEN; i++)
		{
			for (int j = 0; j < NUM_LEVEL; j++)
			{
				std::cout << L[0][j][i] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << "-----------------------" << std::endl;
		for (int i = 0; i < CODE_LEN; i++)
		{
			for (int j = 0; j < NUM_LEVEL; j++)
			{
				std::cout << L[BCT_LAYER][j][i] << " ";
			}
			std::cout << std::endl;
		}
		system("pause");*/


	int upper, lower;

	for (int k = 0; k < Iteration; k++)
	{
		for (int i = 0; i < BCT_LAYER; i++)
		{

			int G = pow(2, i);
			for (int g = 0; g < G; g++)
			{
				int p = pow(2, i + 1);

				for (int j = 0; j < CODE_LEN / p; j++)
				{
					upper = g * (CODE_LEN / G) + j;
					lower = g * (CODE_LEN / G) + (CODE_LEN / p) + j;

					//L[i + 1][upper] = min_sum(L[i][lower]+R[i + 1][lower], L[i][upper]);
					//L[i + 1][lower] = L[i][lower] + min_sum(L[i][upper], R[i + 1][upper]);

					//----update L upper branch
					std::vector<double> temp(NUM_LEVEL);
					double sum = 0;

					for (int l = 0; l < NUM_LEVEL; l++)
					{
						temp[l] = L[i][l][lower] * R[i + 1][l][lower];
						sum += temp[l];
					}
					//-----normal
					for (int l = 0; l < NUM_LEVEL; l++)
					{
						temp[l] /= sum;
					}

					if (NUM_LEVEL == 4)
					{
						L[i + 1][0][upper] = temp[0] * L[i][0][upper] + temp[1] * L[i][1][upper] + temp[2] * L[i][2][upper] + temp[3] * L[i][3][upper];
						L[i + 1][1][upper] = temp[0] * L[i][1][upper] + temp[1] * L[i][0][upper] + temp[2] * L[i][3][upper] + temp[3] * L[i][2][upper];
						L[i + 1][2][upper] = temp[0] * L[i][2][upper] + temp[1] * L[i][3][upper] + temp[2] * L[i][0][upper] + temp[3] * L[i][1][upper];
						L[i + 1][3][upper] = temp[0] * L[i][3][upper] + temp[1] * L[i][2][upper] + temp[2] * L[i][1][upper] + temp[3] * L[i][0][upper];
					}
					else if (NUM_LEVEL == 8)
					{
						L[i + 1][0][upper] = temp[0] * L[i][0][upper] + temp[1] * L[i][1][upper] + temp[2] * L[i][2][upper] + temp[3] * L[i][3][upper] + temp[4] * L[i][4][upper] + temp[5] * L[i][5][upper] + temp[6] * L[i][6][upper] + temp[7] * L[i][7][upper];
						L[i + 1][1][upper] = temp[0] * L[i][1][upper] + temp[1] * L[i][0][upper] + temp[2] * L[i][3][upper] + temp[3] * L[i][2][upper] + temp[4] * L[i][5][upper] + temp[5] * L[i][4][upper] + temp[6] * L[i][7][upper] + temp[7] * L[i][6][upper];
						L[i + 1][2][upper] = temp[0] * L[i][2][upper] + temp[1] * L[i][3][upper] + temp[2] * L[i][0][upper] + temp[3] * L[i][1][upper] + temp[4] * L[i][6][upper] + temp[5] * L[i][7][upper] + temp[6] * L[i][4][upper] + temp[7] * L[i][5][upper];
						L[i + 1][3][upper] = temp[0] * L[i][3][upper] + temp[1] * L[i][2][upper] + temp[2] * L[i][1][upper] + temp[3] * L[i][0][upper] + temp[4] * L[i][7][upper] + temp[5] * L[i][6][upper] + temp[6] * L[i][5][upper] + temp[7] * L[i][4][upper];
						L[i + 1][4][upper] = temp[0] * L[i][4][upper] + temp[1] * L[i][5][upper] + temp[2] * L[i][6][upper] + temp[3] * L[i][7][upper] + temp[4] * L[i][0][upper] + temp[5] * L[i][1][upper] + temp[6] * L[i][2][upper] + temp[7] * L[i][3][upper];
						L[i + 1][5][upper] = temp[0] * L[i][5][upper] + temp[1] * L[i][4][upper] + temp[2] * L[i][7][upper] + temp[3] * L[i][6][upper] + temp[4] * L[i][1][upper] + temp[5] * L[i][0][upper] + temp[6] * L[i][3][upper] + temp[7] * L[i][2][upper];
						L[i + 1][6][upper] = temp[0] * L[i][6][upper] + temp[1] * L[i][7][upper] + temp[2] * L[i][4][upper] + temp[3] * L[i][5][upper] + temp[4] * L[i][2][upper] + temp[5] * L[i][3][upper] + temp[6] * L[i][0][upper] + temp[7] * L[i][1][upper];
						L[i + 1][7][upper] = temp[0] * L[i][7][upper] + temp[1] * L[i][6][upper] + temp[2] * L[i][5][upper] + temp[3] * L[i][4][upper] + temp[4] * L[i][3][upper] + temp[5] * L[i][2][upper] + temp[6] * L[i][1][upper] + temp[7] * L[i][0][upper];
					}
					else
					{
						std::cout << "segment fault";
						system("pause");
					}

					//----limit bound
					for (int j = 0; j < NUM_LEVEL; j++)
					{
						if (L[i + 1][j][upper] < NUMERIC_LIMIT)
							L[i + 1][j][upper] = NUMERIC_LIMIT;
					}


					//---- update L lower branch

					if (NUM_LEVEL == 4)
					{
						L[i + 1][0][lower] = R[i + 1][0][upper] * L[i][0][upper] + R[i + 1][1][upper] * L[i][1][upper] + R[i + 1][2][upper] * L[i][2][upper] + R[i + 1][3][upper] * L[i][3][upper];
						L[i + 1][1][lower] = R[i + 1][0][upper] * L[i][1][upper] + R[i + 1][1][upper] * L[i][0][upper] + R[i + 1][2][upper] * L[i][3][upper] + R[i + 1][3][upper] * L[i][2][upper];
						L[i + 1][2][lower] = R[i + 1][0][upper] * L[i][2][upper] + R[i + 1][1][upper] * L[i][3][upper] + R[i + 1][2][upper] * L[i][0][upper] + R[i + 1][3][upper] * L[i][1][upper];
						L[i + 1][3][lower] = R[i + 1][0][upper] * L[i][3][upper] + R[i + 1][1][upper] * L[i][2][upper] + R[i + 1][2][upper] * L[i][1][upper] + R[i + 1][3][upper] * L[i][0][upper];
					}
					else if (NUM_LEVEL == 8)
					{
						L[i + 1][0][lower] = R[i + 1][0][upper] * L[i][0][upper] + R[i + 1][1][upper] * L[i][1][upper] + R[i + 1][2][upper] * L[i][2][upper] + R[i + 1][3][upper] * L[i][3][upper] + R[i + 1][4][upper] * L[i][4][upper] + R[i + 1][5][upper] * L[i][5][upper] + R[i + 1][6][upper] * L[i][6][upper] + R[i + 1][7][upper] * L[i][7][upper];
						L[i + 1][1][lower] = R[i + 1][0][upper] * L[i][1][upper] + R[i + 1][1][upper] * L[i][0][upper] + R[i + 1][2][upper] * L[i][3][upper] + R[i + 1][3][upper] * L[i][2][upper] + R[i + 1][4][upper] * L[i][5][upper] + R[i + 1][5][upper] * L[i][4][upper] + R[i + 1][6][upper] * L[i][7][upper] + R[i + 1][7][upper] * L[i][6][upper];
						L[i + 1][2][lower] = R[i + 1][0][upper] * L[i][2][upper] + R[i + 1][1][upper] * L[i][3][upper] + R[i + 1][2][upper] * L[i][0][upper] + R[i + 1][3][upper] * L[i][1][upper] + R[i + 1][4][upper] * L[i][6][upper] + R[i + 1][5][upper] * L[i][7][upper] + R[i + 1][6][upper] * L[i][4][upper] + R[i + 1][7][upper] * L[i][5][upper];
						L[i + 1][3][lower] = R[i + 1][0][upper] * L[i][3][upper] + R[i + 1][1][upper] * L[i][2][upper] + R[i + 1][2][upper] * L[i][1][upper] + R[i + 1][3][upper] * L[i][0][upper] + R[i + 1][4][upper] * L[i][7][upper] + R[i + 1][5][upper] * L[i][6][upper] + R[i + 1][6][upper] * L[i][5][upper] + R[i + 1][7][upper] * L[i][4][upper];
						L[i + 1][4][lower] = R[i + 1][0][upper] * L[i][4][upper] + R[i + 1][1][upper] * L[i][5][upper] + R[i + 1][2][upper] * L[i][6][upper] + R[i + 1][3][upper] * L[i][7][upper] + R[i + 1][4][upper] * L[i][0][upper] + R[i + 1][5][upper] * L[i][1][upper] + R[i + 1][6][upper] * L[i][2][upper] + R[i + 1][7][upper] * L[i][3][upper];
						L[i + 1][5][lower] = R[i + 1][0][upper] * L[i][5][upper] + R[i + 1][1][upper] * L[i][4][upper] + R[i + 1][2][upper] * L[i][7][upper] + R[i + 1][3][upper] * L[i][6][upper] + R[i + 1][4][upper] * L[i][1][upper] + R[i + 1][5][upper] * L[i][0][upper] + R[i + 1][6][upper] * L[i][3][upper] + R[i + 1][7][upper] * L[i][2][upper];
						L[i + 1][6][lower] = R[i + 1][0][upper] * L[i][6][upper] + R[i + 1][1][upper] * L[i][7][upper] + R[i + 1][2][upper] * L[i][4][upper] + R[i + 1][3][upper] * L[i][5][upper] + R[i + 1][4][upper] * L[i][2][upper] + R[i + 1][5][upper] * L[i][3][upper] + R[i + 1][6][upper] * L[i][0][upper] + R[i + 1][7][upper] * L[i][1][upper];
						L[i + 1][7][lower] = R[i + 1][0][upper] * L[i][7][upper] + R[i + 1][1][upper] * L[i][6][upper] + R[i + 1][2][upper] * L[i][5][upper] + R[i + 1][3][upper] * L[i][4][upper] + R[i + 1][4][upper] * L[i][3][upper] + R[i + 1][5][upper] * L[i][2][upper] + R[i + 1][6][upper] * L[i][1][upper] + R[i + 1][7][upper] * L[i][0][upper];
					}
					else
					{
						std::cout << "segment fault";
					}

					sum = 0;
					for (int l = 0; l < NUM_LEVEL; l++)
					{
						L[i + 1][l][lower] *= L[i][l][lower];
						sum += L[i + 1][l][lower];
					}
					//----normal
					for (int l = 0; l < NUM_LEVEL; l++)
					{
						L[i + 1][l][lower] /= sum;
					}
					//----limit bound
					for (int j = 0; j < NUM_LEVEL; j++)
					{
						if (L[i + 1][j][lower] < NUMERIC_LIMIT)
							L[i + 1][j][lower] = NUMERIC_LIMIT;
					}

				}
			}

		}

		/*	for (int i = 0; i < CODE_LEN; i++)
			{
				for (int j = 0; j < NUM_LEVEL; j++)
				{
					std::cout << L[0][j][i] << " ";
				}
				std::cout << std::endl;
			}
			std::cout << "-----------------------" << std::endl;
			for (int i = 0; i < CODE_LEN; i++)
			{
				for (int j = 0; j < NUM_LEVEL; j++)
				{
					std::cout << L[BCT_LAYER][j][i] << " ";
				}
				std::cout << std::endl;
			}
			system("pause");*/


		for (int j = 0; j < BCT_LAYER; j++)
		{
			for (int i = BCT_LAYER; i >= 0; i--)
			{

				int G = pow(2, i);
				for (int g = 0; g < G; g++)
				{
					int p = pow(2, i + 1);

					for (int j = 0; j < CODE_LEN / p; j++)
					{
						upper = g * (CODE_LEN / G) + j;
						lower = g * (CODE_LEN / G) + (CODE_LEN / p) + j;

						//R[i][upper] = min_sum(R[i + 1][lower] + L[i][lower], R[i + 1][upper]);
						//R[i][lower] = min_sum(R[i + 1][upper], L[i][upper]) + R[i + 1][lower];

						///////////////////////////////////////////////
						std::vector<double> temp(NUM_LEVEL);
						double sum = 0;

						for (int l = 0; l < NUM_LEVEL; l++)
						{
							temp[l] = R[i + 1][l][lower] * L[i][l][lower];
							sum += temp[l];
						}
						//normal
						for (int l = 0; l < NUM_LEVEL; l++)
						{
							temp[l] /= sum;
						}

						if (NUM_LEVEL == 4)
						{
							R[i][0][upper] = temp[0] * R[i + 1][0][upper] + temp[1] * R[i + 1][1][upper] + temp[2] * R[i + 1][2][upper] + temp[3] * R[i + 1][3][upper];
							R[i][1][upper] = temp[0] * R[i + 1][1][upper] + temp[1] * R[i + 1][0][upper] + temp[2] * R[i + 1][3][upper] + temp[3] * R[i + 1][2][upper];
							R[i][2][upper] = temp[0] * R[i + 1][2][upper] + temp[1] * R[i + 1][3][upper] + temp[2] * R[i + 1][0][upper] + temp[3] * R[i + 1][1][upper];
							R[i][3][upper] = temp[0] * R[i + 1][3][upper] + temp[1] * R[i + 1][2][upper] + temp[2] * R[i + 1][1][upper] + temp[3] * R[i + 1][0][upper];
						}
						else if (NUM_LEVEL == 8)
						{
							R[i][0][upper] = temp[0] * R[i + 1][0][upper] + temp[1] * R[i + 1][1][upper] + temp[2] * R[i + 1][2][upper] + temp[3] * R[i + 1][3][upper] + temp[4] * R[i + 1][4][upper] + temp[5] * R[i + 1][5][upper] + temp[6] * R[i + 1][6][upper] + temp[7] * R[i + 1][7][upper];
							R[i][1][upper] = temp[0] * R[i + 1][1][upper] + temp[1] * R[i + 1][0][upper] + temp[2] * R[i + 1][3][upper] + temp[3] * R[i + 1][2][upper] + temp[4] * R[i + 1][5][upper] + temp[5] * R[i + 1][4][upper] + temp[6] * R[i + 1][7][upper] + temp[7] * R[i + 1][6][upper];
							R[i][2][upper] = temp[0] * R[i + 1][2][upper] + temp[1] * R[i + 1][3][upper] + temp[2] * R[i + 1][0][upper] + temp[3] * R[i + 1][1][upper] + temp[4] * R[i + 1][6][upper] + temp[5] * R[i + 1][7][upper] + temp[6] * R[i + 1][4][upper] + temp[7] * R[i + 1][5][upper];
							R[i][3][upper] = temp[0] * R[i + 1][3][upper] + temp[1] * R[i + 1][2][upper] + temp[2] * R[i + 1][1][upper] + temp[3] * R[i + 1][0][upper] + temp[4] * R[i + 1][7][upper] + temp[5] * R[i + 1][6][upper] + temp[6] * R[i + 1][5][upper] + temp[7] * R[i + 1][4][upper];
							R[i][4][upper] = temp[0] * R[i + 1][4][upper] + temp[1] * R[i + 1][5][upper] + temp[2] * R[i + 1][6][upper] + temp[3] * R[i + 1][7][upper] + temp[4] * R[i + 1][0][upper] + temp[5] * R[i + 1][1][upper] + temp[6] * R[i + 1][2][upper] + temp[7] * R[i + 1][3][upper];
							R[i][5][upper] = temp[0] * R[i + 1][5][upper] + temp[1] * R[i + 1][4][upper] + temp[2] * R[i + 1][7][upper] + temp[3] * R[i + 1][6][upper] + temp[4] * R[i + 1][1][upper] + temp[5] * R[i + 1][0][upper] + temp[6] * R[i + 1][3][upper] + temp[7] * R[i + 1][2][upper];
							R[i][6][upper] = temp[0] * R[i + 1][6][upper] + temp[1] * R[i + 1][7][upper] + temp[2] * R[i + 1][4][upper] + temp[3] * R[i + 1][5][upper] + temp[4] * R[i + 1][2][upper] + temp[5] * R[i + 1][3][upper] + temp[6] * R[i + 1][0][upper] + temp[7] * R[i + 1][1][upper];
							R[i][7][upper] = temp[0] * R[i + 1][7][upper] + temp[1] * R[i + 1][6][upper] + temp[2] * R[i + 1][5][upper] + temp[3] * R[i + 1][4][upper] + temp[4] * R[i + 1][3][upper] + temp[5] * R[i + 1][2][upper] + temp[6] * R[i + 1][1][upper] + temp[7] * R[i + 1][0][upper];
						}
						else
						{
							std::cout << "segment fault";
						}

						//limit bound
						for (int j = 0; j < NUM_LEVEL; j++)
						{
							if (R[i][j][upper] < NUMERIC_LIMIT)
								R[i][j][upper] = NUMERIC_LIMIT;
						}

						/////////////////////////////////////////////

						if (NUM_LEVEL == 4)
						{
							R[i][0][lower] = R[i + 1][0][upper] * L[i][0][upper] + R[i + 1][1][upper] * L[i][1][upper] + R[i + 1][2][upper] * L[i][2][upper] + R[i + 1][3][upper] * L[i][3][upper];
							R[i][1][lower] = R[i + 1][0][upper] * L[i][1][upper] + R[i + 1][1][upper] * L[i][0][upper] + R[i + 1][2][upper] * L[i][3][upper] + R[i + 1][3][upper] * L[i][2][upper];
							R[i][2][lower] = R[i + 1][0][upper] * L[i][2][upper] + R[i + 1][1][upper] * L[i][3][upper] + R[i + 1][2][upper] * L[i][0][upper] + R[i + 1][3][upper] * L[i][1][upper];
							R[i][3][lower] = R[i + 1][0][upper] * L[i][3][upper] + R[i + 1][1][upper] * L[i][2][upper] + R[i + 1][2][upper] * L[i][1][upper] + R[i + 1][3][upper] * L[i][0][upper];
						}
						else if (NUM_LEVEL == 8)
						{
							R[i][0][lower] = L[i][0][upper] * R[i + 1][0][upper] + L[i][1][upper] * R[i + 1][1][upper] + L[i][2][upper] * R[i + 1][2][upper] + L[i][3][upper] * R[i + 1][3][upper] + L[i][4][upper] * R[i + 1][4][upper] + L[i][5][upper] * R[i + 1][5][upper] + L[i][6][upper] * R[i + 1][6][upper] + L[i][7][upper] * R[i + 1][7][upper];
							R[i][1][lower] = L[i][0][upper] * R[i + 1][1][upper] + L[i][1][upper] * R[i + 1][0][upper] + L[i][2][upper] * R[i + 1][3][upper] + L[i][3][upper] * R[i + 1][2][upper] + L[i][4][upper] * R[i + 1][5][upper] + L[i][5][upper] * R[i + 1][4][upper] + L[i][6][upper] * R[i + 1][7][upper] + L[i][7][upper] * R[i + 1][6][upper];
							R[i][2][lower] = L[i][0][upper] * R[i + 1][2][upper] + L[i][1][upper] * R[i + 1][3][upper] + L[i][2][upper] * R[i + 1][0][upper] + L[i][3][upper] * R[i + 1][1][upper] + L[i][4][upper] * R[i + 1][6][upper] + L[i][5][upper] * R[i + 1][7][upper] + L[i][6][upper] * R[i + 1][4][upper] + L[i][7][upper] * R[i + 1][5][upper];
							R[i][3][lower] = L[i][0][upper] * R[i + 1][3][upper] + L[i][1][upper] * R[i + 1][2][upper] + L[i][2][upper] * R[i + 1][1][upper] + L[i][3][upper] * R[i + 1][0][upper] + L[i][4][upper] * R[i + 1][7][upper] + L[i][5][upper] * R[i + 1][6][upper] + L[i][6][upper] * R[i + 1][5][upper] + L[i][7][upper] * R[i + 1][4][upper];
							R[i][4][lower] = L[i][0][upper] * R[i + 1][4][upper] + L[i][1][upper] * R[i + 1][5][upper] + L[i][2][upper] * R[i + 1][6][upper] + L[i][3][upper] * R[i + 1][7][upper] + L[i][4][upper] * R[i + 1][0][upper] + L[i][5][upper] * R[i + 1][1][upper] + L[i][6][upper] * R[i + 1][2][upper] + L[i][7][upper] * R[i + 1][3][upper];
							R[i][5][lower] = L[i][0][upper] * R[i + 1][5][upper] + L[i][1][upper] * R[i + 1][4][upper] + L[i][2][upper] * R[i + 1][7][upper] + L[i][3][upper] * R[i + 1][6][upper] + L[i][4][upper] * R[i + 1][1][upper] + L[i][5][upper] * R[i + 1][0][upper] + L[i][6][upper] * R[i + 1][3][upper] + L[i][7][upper] * R[i + 1][2][upper];
							R[i][6][lower] = L[i][0][upper] * R[i + 1][6][upper] + L[i][1][upper] * R[i + 1][7][upper] + L[i][2][upper] * R[i + 1][4][upper] + L[i][3][upper] * R[i + 1][5][upper] + L[i][4][upper] * R[i + 1][2][upper] + L[i][5][upper] * R[i + 1][3][upper] + L[i][6][upper] * R[i + 1][0][upper] + L[i][7][upper] * R[i + 1][1][upper];
							R[i][7][lower] = L[i][0][upper] * R[i + 1][7][upper] + L[i][1][upper] * R[i + 1][6][upper] + L[i][2][upper] * R[i + 1][5][upper] + L[i][3][upper] * R[i + 1][4][upper] + L[i][4][upper] * R[i + 1][3][upper] + L[i][5][upper] * R[i + 1][2][upper] + L[i][6][upper] * R[i + 1][1][upper] + L[i][7][upper] * R[i + 1][0][upper];
						}
						else
						{
							std::cout << "segment fault";
						}

						sum = 0;
						for (int l = 0; l < NUM_LEVEL; l++)
						{
							R[i][l][lower] *= R[i + 1][l][lower];
							sum += R[i][l][lower];
						}
						//normal
						for (int l = 0; l < NUM_LEVEL; l++)
						{
							R[i][l][lower] /= sum;
						}
						//limit bound
						for (int j = 0; j < NUM_LEVEL; j++)
						{
							if (R[i][j][lower] < NUMERIC_LIMIT)
								R[i][j][lower] = NUMERIC_LIMIT;
						}

						/////////////////////////////////////////////

					}
				}
			}
		}
	}


	/*for (int i = 0; i < CODE_LEN; i++)
	{
		for (int j = 0; j < NUM_LEVEL; j++)
		{
			std::cout << R[BCT_LAYER][j][i] << " ";
		}
		std::cout << std::endl;
	}*/

	//std::cout << "----------";
	std::vector<std::vector<double>> new_app(DATA_LEN, std::vector<double>(NUM_LEVEL));

	for (int i = 0; i < DATA_LEN; i++)
	{
		for (int j = 0; j < NUM_LEVEL; j++)
		{
			new_app[i][j] = L[BCT_LAYER][j][_channel_order_sorted[0][i]];
			//std::cout << new_app[i][j] << " ";
		}
		//std::cout << std::endl;
	}

	for (int i = 0; i < DATA_LEN; i++)
	{
		if (NUM_USER == 1)
		{
			std::cout << "segment fault";
			system("pause");
		}
		else if (NUM_USER == 2)
		{
			//---------- user-A ----------
			if ((new_app[i][0] + new_app[i][1]) <= NUMERIC_LIMIT) appLlr[0][i] = -LLR_LIMIT;
			else if ((new_app[i][2] + new_app[i][3]) <= NUMERIC_LIMIT) appLlr[0][i] = LLR_LIMIT;
			else appLlr[0][i] = log((new_app[i][0] + new_app[i][1]) / (new_app[i][2] + new_app[i][3]));
			//---------- user-B ----------
			if ((new_app[i][0] + new_app[i][2]) <= NUMERIC_LIMIT) appLlr[1][i] = -LLR_LIMIT;
			else if ((new_app[i][1] + new_app[i][3]) <= NUMERIC_LIMIT) appLlr[1][i] = LLR_LIMIT;
			else appLlr[1][i] = log((new_app[i][0] + new_app[i][2]) / (new_app[i][1] + new_app[i][3]));
		}
		else if (NUM_USER == 3)
		{
			//---------- user-A ----------
			if ((new_app[i][0] + new_app[i][1] + new_app[i][2] + new_app[i][3]) <= NUMERIC_LIMIT) appLlr[0][i] = -LLR_LIMIT;
			else if ((new_app[i][4] + new_app[i][5] + new_app[i][6] + new_app[i][7]) <= NUMERIC_LIMIT) appLlr[0][i] = LLR_LIMIT;
			else appLlr[0][i] = log((new_app[i][0] + new_app[i][1] + new_app[i][2] + new_app[i][3]) / (new_app[i][4] + new_app[i][5] + new_app[i][6] + new_app[i][7]));
			//---------- user-B ----------
			if ((new_app[i][0] + new_app[i][1] + new_app[i][4] + new_app[i][5]) <= NUMERIC_LIMIT) appLlr[1][i] = -LLR_LIMIT;
			else if ((new_app[i][2] + new_app[i][3] + new_app[i][6] + new_app[i][7]) <= NUMERIC_LIMIT) appLlr[1][i] = LLR_LIMIT;
			else appLlr[1][i] = log((new_app[i][0] + new_app[i][1] + new_app[i][4] + new_app[i][5]) / (new_app[i][2] + new_app[i][3] + new_app[i][6] + new_app[i][7]));
			//---------- user-C ----------
			if ((new_app[i][0] + new_app[i][2] + new_app[i][4] + new_app[i][6]) <= NUMERIC_LIMIT) appLlr[2][i] = -LLR_LIMIT;
			else if ((new_app[i][1] + new_app[i][3] + new_app[i][5] + new_app[i][7]) <= NUMERIC_LIMIT) appLlr[2][i] = LLR_LIMIT;
			else appLlr[2][i] = log((new_app[i][0] + new_app[i][2] + new_app[i][4] + new_app[i][6]) / (new_app[i][1] + new_app[i][3] + new_app[i][5] + new_app[i][7]));
		}

	}

	//	std::cout << std::endl << "------------";

	for (int i = 0; i < NUM_USER; i++)
	{
		for (int j = 0; j < DATA_LEN; j++)
		{
			decodedResult[i][j] = HARD(appLlr[i][j]);
			//decodedResult[i][j] = appLlr[i][j];
			//std::cout << decodedResult[i][j] << " ";
		}
		//std::cout << std::endl;
	}
	//	system("pause");
}

std::vector<int> PolarCode::decode_BP_sep(double*  llr, int nuser)
{
	for (int i = 0; i < _block_length; i++)
	{
		if (isinf(llr[i]))
		{
			llr[i] = signbit(llr[i]) ? -40 : 40;
		}
		llr[i] > 40 ? llr[i] = 40 : 0;
		llr[i] < -40 ? llr[i] = -40 : 0;
	}
	std::vector<std::vector<double>> R(BCT_LAYER+1, std::vector<double>(CODE_LEN));
	std::vector<std::vector<double>> L(BCT_LAYER+1, std::vector<double>(CODE_LEN));
	for (int i = 0; i < CODE_LEN; i++)
	{
		if (_frozen_bits[nuser][i])
		{
			R[BCT_LAYER][i] = 40;
		}
		L[0][i] = llr[i];
	}
	//system("pause");

	int upper, lower;


	for (int k = 0; k < Iteration; k++)
	{
		for (int i = 0; i < BCT_LAYER; i++)
		{

			int G = pow(2, i);
			for (int g = 0; g < G; g++)
			{
				int p = pow(2, i + 1);

				for (int j = 0; j < CODE_LEN / p; j++)
				{
					upper = g * (CODE_LEN / G) + j;
					lower = g * (CODE_LEN / G) + (CODE_LEN / p) + j;

					L[i + 1][upper] = min_sum(L[i][lower] + R[i + 1][lower], L[i][upper]);
					L[i + 1][lower] = L[i][lower] + min_sum(L[i][upper], R[i + 1][upper]);

					L[i + 1][upper] = L[i + 1][upper] < -LLR_LIMIT ? -LLR_LIMIT : L[i + 1][upper] > LLR_LIMIT ? LLR_LIMIT : L[i + 1][upper];
					L[i + 1][lower] = L[i + 1][upper] < -LLR_LIMIT ? -LLR_LIMIT : L[i + 1][lower] > LLR_LIMIT ? LLR_LIMIT : L[i + 1][lower];
				}
			}

		}
		
		for (int j = 0; j < BCT_LAYER; j++)
		{
			for (int i = BCT_LAYER; i >= 0; i--)
			{

				int G = pow(2, i);
				for (int g = 0; g < G; g++)
				{
					int p = pow(2, i + 1);

					for (int j = 0; j < CODE_LEN / p; j++)
					{
						upper = g * (CODE_LEN / G) + j;
						lower = g * (CODE_LEN / G) + (CODE_LEN / p) + j;

						R[i][upper] = min_sum(R[i + 1][lower] + L[i][lower], R[i + 1][upper]);
						R[i][lower] = min_sum(R[i + 1][upper], L[i][upper]) + R[i + 1][lower];

						R[i + 1][upper] = R[i + 1][upper] < -LLR_LIMIT ? -LLR_LIMIT : R[i + 1][upper] > LLR_LIMIT ? LLR_LIMIT : R[i + 1][upper];
						R[i + 1][lower] = R[i + 1][upper] < -LLR_LIMIT ? -LLR_LIMIT : R[i + 1][lower] > LLR_LIMIT ? LLR_LIMIT : R[i + 1][lower];
					}
				}
			}
		}
	}

	std::vector<int> message_output;

	for (int i = 0; i < CODE_LEN; i++)
	{
		if (_frozen_bits[nuser][i])
			continue;
		message_output.push_back(HARD(L[BCT_LAYER][i]));
	}
	
	return message_output;
}

double PolarCode::min_sum(double a, double b)
{
	/*if (a * b >= 0)
		return std::min(std::abs(a), std::abs(b));
	else
		return -(std::min(std::abs(a), std::abs(b)));*/

	if (40 > std::max(std::abs(a), std::abs(b)))
		return std::log((exp(a + b) + 1) / (exp(a) + exp(b)));
	else
		return  (double)((a < 0) ? -1 : (a > 0)) * ((b < 0) ? -1 : (b > 0)) * std::min(std::abs(a), std::abs(b));
}

std::vector<double> PolarCode::vector_product(std::vector<double> v1, std::vector<double> v2)
{
	std::vector<double> re;
	std::vector<double> temp(NUM_LEVEL);
	double sum = 0;

	for (int i = 0; i < NUM_LEVEL; i++)
	{
		temp[i] =v1[i] * v2[i];
		sum += temp[i];
	}
	//normal
	for (int l = 0; l < NUM_LEVEL; l++)
	{
		temp[l] /= sum;
		re.push_back(temp[l]);
	}

	return re;

}


