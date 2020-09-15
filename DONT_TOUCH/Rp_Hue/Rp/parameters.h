#pragma once

#include "ldpc.h"
#include "lt.h"

//---------- Raptor code setting ----------

#define		LDPC_H_COL			2000
#define		LDPC_H_ROW			100
#define		INFILENAME			"H_2000_100.txt"

//#define		LDPC_H_COL			4000
//#define		LDPC_H_ROW			200
//#define		INFILENAME			"H_4000_200.txt"

//#define		LDPC_H_COL			10000
//#define		LDPC_H_ROW			500
//#define		INFILENAME			"H_10000_500.txt"

//#define		LDPC_H_COL			10000
//#define		LDPC_H_ROW			200
//#define		INFILENAME			"H_10000_200.txt"

//#define		LDPC_H_COL			12000
//#define		LDPC_H_ROW			600
//#define		INFILENAME			"H_12000_600.txt"

#define		DEC_SCHEME			1														// 1: tandem decoding; 2: joint decoding
#define		LDPC_IT				10														// number of iterations for LDPC code
#define		LT_IT				20	 													// number of iterations for LT code
#define		JOINT_IT			30														// number of iterations for joint decoding

#define		IR_LEN				40														// size of retransmissions
#define		IR_MAX				100														// maximum number of retransmissions

#define		DEGREE_SET			312														// output check node degree distribution
#define		SYSTEMATIC			1														// 1: systematic; 0: non-systematic

#define		DATA_LEN			( LDPC_H_COL - LDPC_H_ROW )								// data length
#define		INTER_CODE_LEN		LDPC_H_COL												// length of LDPC-coded sequence
#define		MAX_CODE_LEN		( INTER_CODE_LEN + ( IR_LEN * IR_MAX ) )				// maximun length of LT-coded sequence

//---------- simulator setting ----------

#define		BLOCK_NUM			100													// number of blocks to be simulated 

#define		SNR_NUM				5														// number of SNR points to be simulated
#define		SNR_START			1														// in dB
#define		SNR_STEP			1														// in dB

#define		HARD(x)				( (x) > 0 ? 0 : 1 )			

//---------- function prototype ----------

void	Channel(double stdDev, int *tx, double *rx);
double	Capacity(double variance);

void	RaptorEncoder(LDPC &ldpc, LT &lt, int *data, int *inter_code, int *tx);
void	TandemDecoder(LT &lt, LDPC &ldpc, double variance, double *rx, double *inter_llr, double *decoded_llr, int *data, int &ack, int &ir_count);
void	JointDecoder(LT &lt, LDPC &ldpc, double variance, double *rx, double *decoded_llr, int *data, int &ack, int &ir_count);
void	LLR_C2V(LDPC &ldpc, LT &lt, int current_len);
void	LLR_V2C(LDPC &ldpc, LT &lt, int current_len);