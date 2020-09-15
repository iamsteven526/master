#pragma once
#include <vector>
#include "lt.h"
#include "ldpc.h"
using namespace std;

class Raptor
{
public:

	Raptor(LDPC&,LT&,int,int,int);
//	vector<int> encoder(vector<int>);
	vector<double> joint_decoder(int, double, vector<double> , bool&,LDPC&);
	vector<double> tandem_decoder(double, vector<double>, int, bool&, LDPC&, LT&);
	vector<int> row_address;
	vector<int> column_address;
	vector<int> column_bits;
	vector<int> row_bits;

private:
	int data_len;
	int inter_code_len;
	int max_code_len;

};
