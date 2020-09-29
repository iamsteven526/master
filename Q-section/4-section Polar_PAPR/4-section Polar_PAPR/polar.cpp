

#include <iostream>
#include <random>
#include <algorithm>
#include <functional>
#include <fstream>
#include "polar.h"
#include "parameters.h"

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
			system("pause");
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
		std::cout << i << ". ";
		for (int j = 0; j < _block_length; j++)
		{
			std::cout << G[i][j] << " ";
		}
		std::cout << std::endl;
	}
	system("pause");*/
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



bool PolarCode::crc_check(uint8_t *info_bit_padded, int nuser)
{
	bool crc_pass = true;
	for (uint16_t i = _info_length; i < _info_length + _crc_size; ++i)
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
	}
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



std::vector<int> PolarCode::decode_scl(int nuser)
{
	for (uint16_t phi = 0; phi < _block_length; ++phi)
	{
		if (_llr_based_computation)
		{
			recursivelyCalcLLR(_n, phi);  //修改message passing的方式
		}
		else
		{
			recursivelyCalcP(_n, phi);  //目前用不到
		}
		if (_frozen_bits.at(nuser).at(phi) == 1)
		{
			continuePaths_FrozenBit(phi, 0);
		}
		else
		{
			continuePaths_UnfrozenBit(phi, 0);
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
		_arrayPointer_Info.at(s) = new uint8_t[_block_length]();
		for (uint16_t lambda = 0; lambda < _n + 1; ++lambda)
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
	//std::cout << "ok"; system("pause");
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
	if (_llr_based_computation)
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



void PolarCode::recursivelyCalcP(uint16_t lambda, uint16_t phi)
{
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
			sigma = (std::max)(sigma, p_lambda[2 * beta]);
			sigma = (std::max)(sigma, p_lambda[2 * beta + 1]);
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
		double* llr_lambda = getArrayPointer_LLR(lambda, l);
		double* llr_lambda_1 = getArrayPointer_LLR(lambda - 1, l);
		uint8_t* c_lambda = getArrayPointer_C(lambda, l);
		for (uint16_t beta = 0; beta < (1 << (_n - lambda)); ++beta)
		{
			if ((phi % 2) == 0) // upper branch
			{					// min sum
				if (40 > (std::max)(std::abs(llr_lambda_1[2 * beta]), std::abs(llr_lambda_1[2 * beta + 1])))
				{
					llr_lambda[beta] = std::log((exp(llr_lambda_1[2 * beta] + llr_lambda_1[2 * beta + 1]) + 1) /
						(exp(llr_lambda_1[2 * beta]) + exp(llr_lambda_1[2 * beta + 1])));
				}
				else
				{
					llr_lambda[beta] = (double)((llr_lambda_1[2 * beta] < 0) ? -1 : (llr_lambda_1[2 * beta] > 0)) *
						((llr_lambda_1[2 * beta + 1] < 0) ? -1 : (llr_lambda_1[2 * beta + 1] > 0)) *
						(std::min)(std::abs(llr_lambda_1[2 * beta]), std::abs(llr_lambda_1[2 * beta + 1]));
				}
			}
			else // lower branch
			{
				uint8_t u_p = c_lambda[2 * beta];
				llr_lambda[beta] = (1 - 2 * u_p) * llr_lambda_1[2 * beta] + llr_lambda_1[2 * beta + 1];
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
		uint8_t* c_lambda = getArrayPointer_C(lambda, l);
		uint8_t* c_lambda_1 = getArrayPointer_C(lambda - 1, l);
		for (uint16_t beta = 0; beta < (1 << (_n - lambda)); ++beta)
		{
			c_lambda_1[2 * (2 * beta) + (psi % 2)] = (uint8_t)((c_lambda[2 * beta] + c_lambda[2 * beta + 1]) % 2);
			c_lambda_1[2 * (2 * beta + 1) + (psi % 2)] = c_lambda[2 * beta + 1];
		}
	}
	if ((psi % 2) == 1)
	{
		recursivelyUpdateC((uint16_t)(lambda - 1), psi);
	}
}

void PolarCode::continuePaths_FrozenBit(uint16_t phi, int nuser)
{
	for (uint16_t l = 0; l < _list_size; ++l)
	{
		if (_activePath.at(l) == 0)
		{
			continue;
		}

		uint8_t* c_m = getArrayPointer_C(_n, l);
		//---- frozen value assumed to be zero
		c_m[(phi % 2)] = 0;
		if (_llr_based_computation)
		{
			double* llr_p = getArrayPointer_LLR(_n, l);
			_pathMetric_LLR.at(l) += log(1 + exp(-llr_p[0]));
		}
		_arrayPointer_Info.at(l)[phi] = 0;
	}
}

void PolarCode::continuePaths_UnfrozenBit(uint16_t phi, int nuser)
{
	std::vector<double> probForks((unsigned long)(2 * _list_size));
	std::vector<double> probabilities;
	std::vector<uint8_t> contForks((unsigned long)(2 * _list_size));
	uint16_t i = 0;
	for (unsigned l = 0; l < _list_size; ++l)
	{
		if (_activePath.at(l) == 0)
		{
			probForks.at(2 * l) = NAN;
			probForks.at(2 * l + 1) = NAN;
		}
		else
		{
			if (_llr_based_computation)
			{
				double* llr_p = getArrayPointer_LLR(_n, l);
				probForks.at(2 * l) = -(_pathMetric_LLR.at(l) + log(1 + exp(-llr_p[0])));
				probForks.at(2 * l + 1) = -(_pathMetric_LLR.at(l) + log(1 + exp(llr_p[0])));
			}
			else
			{
				double *p_m = getArrayPointer_P(_n, l);
				probForks.at(2 * l) = p_m[0];
				probForks.at(2 * l + 1) = p_m[1];
			}
			probabilities.push_back(probForks.at(2 * l));
			probabilities.push_back(probForks.at(2 * l + 1));
			i++;
		}
	}
	uint16_t rho = _list_size;
	if ((2 * i) < _list_size)
	{
		rho = (uint16_t)2 * i;
	}
	for (uint8_t l = 0; l < 2 * _list_size; ++l)
	{
		contForks.at(l) = 0;
	}
	std::sort(probabilities.begin(), probabilities.end(), std::greater<double>());
	double threshold = probabilities.at((unsigned long)(rho - 1));
	uint16_t num_paths_continued = 0;
	for (uint8_t l = 0; l < 2 * _list_size; ++l)
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
		for (uint8_t l = 0; l < 2 * _list_size; ++l)
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
		if (contForks.at(2 * l) == 0 && contForks.at(2 * l + 1) == 0)
		{
			killPath(l);
		}
	}
	for (unsigned l = 0; l < _list_size; ++l)
	{
		if (contForks.at(2 * l) == 0 && contForks.at(2 * l + 1) == 0)
		{
			continue;
		}
		uint8_t* c_m = getArrayPointer_C(_n, l);
		if (contForks.at(2 * l) == 1 && contForks.at(2 * l + 1) == 1)
		{
			c_m[(phi % 2)] = 0;
			uint16_t l_p = clonePath(l);
			c_m = getArrayPointer_C(_n, l_p);
			c_m[(phi % 2)] = 1;
			std::copy(_arrayPointer_Info.at(l), _arrayPointer_Info.at(l) + phi, _arrayPointer_Info.at(l_p));
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
		else
		{
			if (contForks.at(2 * l) == 1)
			{
				c_m[(phi % 2)] = 0;
				_arrayPointer_Info.at(l)[phi] = 0;
				if (_llr_based_computation)
				{
					double* llr_p = getArrayPointer_LLR(_n, l);
					_pathMetric_LLR.at(l) += log(1 + exp(-llr_p[0]));
				}
			}
			else
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
}

uint16_t PolarCode::findMostProbablePath(bool check_crc, int nuser)
{
	uint16_t l_p = 0;
	double p_p1 = 0;
	double p_llr = (std::numeric_limits<double>::max)();
	bool path_with_crc_pass = false;
	for (uint16_t l = 0; l < _list_size; ++l)
	{
		if (_activePath.at(l) == 0)
		{
			continue;
		}
		if ((check_crc) && (!crc_check(_arrayPointer_Info.at(l), nuser)))
		{
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
		else
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





double PolarCode::min_sum(double a, double b)
{
	if (a * b >= 0)
		return (std::min)(std::abs(a), std::abs(b));
	else
		return -((std::min)(std::abs(a), std::abs(b)));
}




