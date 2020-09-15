#pragma once
#include<vector>
#include"ldpc.h"  
#include"lt.h"    
#include"raptor.h"


using namespace std;
#define		LDPC_H_ROW			100
#define		LDPC_H_COL			2000
#define		INFILENAME			"H_2000_100.txt"


//#define		LDPC_H_ROW			504
//#define		LDPC_H_COL			1008
//#define		INFILENAME			"H_1008_504.txt"

#define		DEC_SCHEME			2														// 1: tandem decoding; 2: joint decoding
#define		LDPC_IT				10														// number of iterations for LDPC code
#define		LT_IT				20	 													// number of iterations for LT code
#define		JOINT_IT			30														// number of iterations for joint decoding

#define		IR_LEN				40														// size of retransmissions
#define		IR_MAX				60														// maximum number of retransmissions

#define		DEGREE_SET			312													// output check node degree distribution
#define		SYSTEMATIC			1														// 1: systematic; 0: non-systematic

#define		DATA_LEN			( LDPC_H_COL - LDPC_H_ROW )							        // data length
#define		INTER_CODE_LEN		LDPC_H_COL												// length of LDPC-coded sequence
#define		MAX_CODE_LEN		( INTER_CODE_LEN + ( IR_LEN * IR_MAX ) )				// maximun length of LT-coded sequence

//---------- simulator setting ----------

#define		BLOCK_NUM			100													// number of blocks to be simulated 

#define		SNR_NUM				5														// number of SNR points to be simulated
#define		SNR_START			1														// in dB
#define		SNR_STEP			1														// in dB

#define		HARD(x)				( (x) > 0 ? 0 : 1 )			

// function prototype
void channel(double, vector<int>&, vector<double>&,int);
