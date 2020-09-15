#pragma once
#include <vector>
 
#define e 2.718281828
using namespace std;
class LDPC 
{
public:

	LDPC(int,int);
	vector<int> encoder(vector<int>);
	vector<double> decoder(double,vector<double> &,bool&);
	vector<int> row_address;
	vector<int> column_address;
	vector<int> column_bits;
	vector<int> row_bits;
	int total_number;

private:
	
	int column;
	int row;

	vector<int> message;
	vector<int> codeword;
	vector<vector<int>> G;
	vector<vector<int>> H;




};