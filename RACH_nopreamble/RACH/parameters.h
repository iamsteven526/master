#pragma once

#include "ldpc.h"

//---------- LDPC code ----------

#define		LDPC_H_COL				1008											// number of columns
#define		LDPC_H_ROW				504												// number of rows

#define		INFILENAME				"H_1008_504.txt"								// (3,6)-regular
//#define		INFILENAME				"H_1008_504_diff_met_dv_12.txt"					// designed for differential encoding

#define		LDPC_IT					150												// maximum number of decoding iterations

#define		DATA_LEN				( LDPC_H_COL - LDPC_H_ROW )						// data length
#define		CODE_LEN				( LDPC_H_COL )									// codeword length

#define		JCD						0												// joint channel decoding; 1: enable, 0: disable

//---------- channel model ----------

#define		CH_TYPE					1												// frequency-selective Rayleigh fading
																					// 1: quasic-static
																					// 2: time-varying
#define		K						0												// 1: K=0 ; Rayleight Channel

#define		TAP_NUM					5												// number of channel taps

#define		CARRIER_FREQ			2												// carrier frequency (GHz)
#define		CHIP_RATE				1.92 //30.72									// sampling frequency (MHz)

#define		UE_SPEED				150												// mobile speed (km/hour)
#define		TEST_FRAME_SIZE			2												// the maximum packets in a frame
#define		SLOTED					0												// 1: enable, 0: disable
#define		ALOHA					1												// 1: Aloha + GDMA,  0: Aloha
#define		PREABLE					0											    // 1: enable, 0: disable	
#define		PREABLE_LEN				63												
#define		Time_Estimate			0												// 0: disable, 1: enable
#define		Phase_Estimate			0												// 0: disable, 1: enable(first channel estimation) change if segfault
#define		Root_Num				1
#define		Threshold				0.15

//----------  pulse shaping ----------

#define		Unit					1												// units in one chip
#define		UP_RATE					Unit											// up-sampling rate
#define		TRUNCATION				10												// truncation of shaping filter		  							
#define		ROLL_OFF				0.22											// roll-off factor
#define		DRIFT_RANGE				6												// drifting range
#define		OVER_SAMPLE			    0												// 1: enable 0: disable
#define		PULSE_SHAPE				0												// 1: enable 0: disable
#define		OVER_SAMPLE_RATE		4


//---------- system ----------

#define		NUM_USER				10												// number of users

#define		NUM_LEVEL				( 1 << NUM_USER )								// number of levels of superimposed signal

#define		FFT_LAYER				4												// number of FFT layers
#define		FFT_POINT				( 1 << FFT_LAYER )								// number of FFT points
#define		FFT_SEGMENT				( CODE_LEN / FFT_POINT )						// codeword is divided into multiple segments for OFDM transmission
#define		CP_LEN					( TAP_NUM - 1 + 4)									// length of cyclic prefix
#define		CS_LEN					0
#define		CP_TYPE					1												// 1: cp, 0:cp + cs

//---------- differential encoding ----------

#define		DIFF_DEC				0												// differential decoding; 1: enable, 0: disable

#define		TURBO_DEC				0												// turbo processing; 1: enable, 0: disable
#define		TURBO_IT				5												// number of turbo iterations

#define		JOINT_DEC				0												// joint decoding; 1: enable, 0: disable
#define		JOINT_IT				100												// number of joint decoding iterations
#define		INNER_IT				1												// number of inner iterations

#define		DIFF_ENC				( DIFF_DEC || TURBO_DEC || JOINT_DEC )			// differential encoding; 1: enable, 0: disable

//---------- channel estimation ----------

#define		CE_SCHEME				1												// 0: cluster-based, 1: ideal

#define		CE_METHOD				1												// 0: preamble estimation only;  1: preamble estimation + modified k-means++

#define		GROUP_SIZE				2												// number of groups to be classified

#define		EM_GMM					1												// EM clustering with Gaussian mixture model; 1: enable, 0: disable

#define		SLIDING					0												// sliding window; 1: enable, 0: disable
#define		WINDOW_SIZE				11												// window size (odd)
#define		EQULIZER				1												// 1: enable 2: disable
#define		MODIFIED_LLR			1												// 1: enable, 0: disable
#define		MODIFIER				0.01		

//---------- simulator setting ----------

#define		NUMERIC_LIMIT			1e-100
#define		LLR_LIMIT				210

#define		BLOCK_NUM				50000											// number of blocks to be simulated 

#define		G_NUM					15												// number of SNR points to be simulated
#define		G_START				    1.1											// in dB
#define		G_STEP					0.2												// in dB

#define		HARD(x)					( (x) > 0 ? 0 : 1 )

//--------------------

void	Encoder(LDPC &ldpc, int ***data, int ***codeword, int *packet_num);
void	DiffEncoding(int *codeword);
void	Modulator(int ***codeword, double ****chip, int* packet_num);
void	MultiCarrierMapper(double ****chip, double *****tx, int *packet_num);
void	MultiCarrierDemapper(double *****ClassifyRx, double *****postRx, int *packet_num);
void	Detector(LDPC &ldpc, int ***data, double ***appLlr, double ***refLlr, long double *errCount, int *packet_num, bool **estimate_packet_time, double **segment_num, double **segment_error, int** packet_time);
void	DiffDecoding(double *appLlr, double *refLlr);

void	UpSampling(double *****tx, double *****preTx, int* packet_num);
double	SquareRootRaisedCosine(double m);
void	SRRCGeneration(double *srrc, double drift);
void	PulseShaping(double *****preTx, double *****tx, double *txFilter, double *drift, int *packet_num, int **packet_time);
void	ReceivingFilter(double **rx, double ***tempRx, double *rxFilter, int *packet_num, int **packet_time);

void	EnergyProfile(double *****h, double *****H, int *packet_num, int **packet_time, double*** h_preamble, double*** H_preamble);
void	MultipleAccessChannel(double stdDev, double *****h, double *****tx, double **rx, double *drift, double *****H, int *packet_num, int **packet_time, double *txFilter, double*** ZC_sequence, double*** h_preamble);

void	MLDT(LDPC &ldpc, double stdDev, double *****H, double *****postRx, double ****app, double ***appLlr, double ***refLlr, double *****estimate, int *packet_num, int***Cluster_num, int ****Cluster_gain, double snrdB);

void	FFT(double *pr, double *pi, int n, int k, double *fr, double *fi, int l, int il);

void	Clustering(double *****postRx, double variance, double *****estimate, double *****H, int ***cluster_num, int *packet_num, double &mse, bool** estimate_packet_time, double** cluster_sample, double **centroid, int **group, int *groupSize, double *variation, int **packet_time, double ** softAssign);
void	Grouping(double** cluster_sample, double** centroid, int** group, int* groupSize, int CLUSTER_USER, int CLUSTER_GROUP, int CLUSTER_LEN);
double	EuclideanDistance(double* x, double* y);
void	CentroidRenewal(double** cluster_sample, double** centroid, int** group, int* groupSize, double* variation, int CLUSTER_GROUP);
bool	ConditionCheck(double** centroid, int* groupSize, int CLUSTER_GROUP);
void	EMClustering(double** rx, double** centroid, double** softAssign, double variance, int CLUSTER_LEN, int CLUSTER_GROUP);

void	Packet_generater(int *packet_num, int **packet_time, long double& packet_sum, double packet_rate, long double* errCount);
void	ReceivingClassify(double **rx, double *****ClassifyRx, int *packet_num, int **packet_time, int*** Cluster_num, int**** Cluster_gain);
void	Time_Num_Estimation(double **rx, bool** estimate_packet_time, double& Time_deviation, double &Time_error,int* packet_num,int** packet_time, int *Time_Offset_Error, double***** estimate, double*** ZC_sequence);
void	ZC_Generater(double ***ZC_sequence);

double  LogMaxFunction(double x, double y);



template <class T>
void Swaping(T **x, T **y)
{
	T *temp = *x;
	*x = *y;
	*y = temp;
}