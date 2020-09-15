#pragma once
#include<vector>
//#include"parameter.h"
using namespace std;
class LT
{
public:
	LT(int,int);
	vector<int> encoder(vector<int>,int);
	vector<double> decoder(double,vector<double>,int);
	vector<int> row_address;
	vector<int> column_address;
	vector<int> row_bits;
	vector<int> column_bits;
	int total_number;
	
private:
	
	int inter_code_len;
	int max_code_len;


	vector<vector<int>> G;
	
	

};