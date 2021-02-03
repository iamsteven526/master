#pragma once

#include "ldpc.h"
#include "polar.h"

#define		CH_CODING_TYPE			0												// 1: LDPC ; 0: Polar 
//---------- Polar code ----------
#define		BCT_LAYER				10												// number of layers in channel polarization

#define		INF						1000

#define		CRC_LEN					4												// number of bits in CRC
#define		LIST_SIZE				32												// list size 

#define		PUNCTURING				0												// 1: enable, 0: disable
#define		PUNCTURED_BIT			16												// number of bits to be punctured

#define		PCC_METHOD				4												// 0: Bhattacharyya bound
																					// 1: Capacity bound
																					// 2: Piecewise linear approximation
																					// 3: Gaussain approximation
																					// 4: 3GPP; fixed lengh of 1024
																					// 5: Section NBC; fixed lengh of 1024

#define		DESIGHED_SNR			1												// SNR per transmitted symbol
#define		POLAR_DECODING_TYPE		1												// 1: SCL,  0: BP
#define		JOINT					0												// 1: JCD,  0: SCD
#define		Iteration				50
#define		NBC						0												// Using the structure of RM code
#define		INTERLEAVER				1												// 1: enable, 0 : disable
#define		INTERLEAVER_type		0												// 1: random, 0 : non-random

#define		DATA_LEN				(512 + CRC_LEN - (CRC_LEN * POLAR_DECODING_TYPE))			// data length  527 for NBC  , 512 for CSI
#define		CODE_LEN				( 1 << BCT_LAYER )								// codeword length
//---------- LDPC code ----------

#define		LDPC_H_COL				1008											// number of columns

#define		LDPC_H_ROW				504												// number of rows

#define		INFILENAME				"H_1008_504.txt"								// (3,6)-regular
//#define		INFILENAME				"H_1008_504_diff_met_dv_12.txt"					// designed for differential encoding
#define		LDPC_IT					100												// maximum number of decoding iterations

#define		JCD						1												// joint channel decoding; 1: enable, 0: disable

//#define		DATA_LEN				( LDPC_H_COL - LDPC_H_ROW )						// data length
//#define		CODE_LEN				( LDPC_H_COL )									// codeword length


	//---------- differential encoding ----------

	#define		DIFF_DEC				0												// differential decoding; 1: enable, 0: disable

	#define		TURBO_DEC				0												// turbo processing; 1: enable, 0: disable
	#define		TURBO_IT				5												// number of turbo iterations

	#define		JOINT_DEC				0												// joint decoding; 1: enable, 0: disable
	#define		JOINT_IT				100												// number of joint decoding iterations
	#define		INNER_IT				1												// number of inner iterations

	#define		DIFF_ENC				( DIFF_DEC || TURBO_DEC || JOINT_DEC )			// differential encoding; 1: enable, 0: disable

//---------- channel model ----------

#define		CH_TYPE					1												// frequency-selective Rayleigh fading
																					// 1: quasic-static
																					// 2: time-varying
#define		K						0												// 1: K=0 ; Rayleight Channel

#define		TAP_NUM					5												// number of channel taps

#define		CARRIER_FREQ			2												// carrier frequency (GHz)
#define		CHIP_RATE				1.92 //30.72									// sampling frequency (MHz)

#define		UE_SPEED				150												// mobile speed (km/hour)

//---------- synchronization ----------

#define		SYNCHRONOUS				1												// 1: synchronous, 0: unknown timing drifts

#define		UP_RATE					4												// up-sampling rate

#define		TRUNCATION				6												// truncation of shaping filter		  							
#define		ROLL_OFF				0.22											// roll-off factor

#define		DRIFT_RANGE				6												// drifting range

#if SYNCHRONOUS == 0
	#define		OVER_SAMPLE			    1												// 1: enable 0: disable
	#define		PULSE_SHAPE				1												// 1: enable 0: disable
	#define		CP_TYPE					0												// 1: cp, 0:cp + cs
#else
	#define		OVER_SAMPLE			    0												// 1: enable 0: disable
	#define		PULSE_SHAPE				0												// 1: enable 0: disable
	#define		CP_TYPE					1												// 1: cp, 0:cp + cs
#endif 

#define		OVER_SAMPLE_RATE		4


//---------- system ----------

#define		NUM_USER				4												// number of users

#define		NUM_LEVEL				( 1 << NUM_USER )								// number of levels of superimposed signal

#define		FFT_LAYER				5												// number of FFT layers
#define		FFT_POINT				( 1 << FFT_LAYER )								// number of FFT points
#define		FFT_SEGMENT				( CODE_LEN / FFT_POINT )						// codeword is divided into multiple segments for OFDM transmission
#define		CP_LEN					( TAP_NUM - 1)									// length of cyclic prefix
#define		CS_LEN					4

//---------- channel estimation ----------

#define		CE_SCHEME				1												// 0: cluster-based, 1: ideal

#define		INI_METHOD				2												// 0: LBG, 1: k-means++, 2: modified k-means++

#define		GROUP_SIZE				NUM_LEVEL										// number of groups to be classified

#define		EM_GMM					1												// EM clustering with Gaussian mixture model; 1: enable, 0: disable

#define		SLIDING					0												// sliding window; 1: enable, 0: disable
#define		WINDOW_SIZE				31												// window size (odd)

//---------- simulator setting ----------

#define		NUMERIC_LIMIT			1e-100
#define		LLR_LIMIT				230

#define		BLOCK_NUM				3000000										// number of blocks to be simulated 

#define		SNR_NUM					9												// number of SNR points to be simulated
#define		SNR_START				16												// in dB
#define		SNR_STEP				2												// in dB

#define		HARD(x)					( (x) > 0 ? 0 : 1 )

//--------------------

void	Encoder(LDPC &ldpc, PolarCode &polar,int **data, int **codeword, int **Interleaver);
void	DiffEncoding(int *codeword);
void	Modulator(int **codeword, double ***chip);
void	MultiCarrierMapper(double ***chip, double ****tx);
void	MultiCarrierDemapper(double ***rx, double ***postRx, double* drift);
void	Detector(LDPC &ldpc, PolarCode &polar, int **data, double **appLlr, double **refLlr, long double *errCount, double **app, int** Interleaver,int *error_bits_count);
void	DiffDecoding(double *appLlr, double *refLlr);

void	TrellisConstruction(int ***trellis);
void	BCJR(int ***trellis, double *refLlr, double *rxLlr, double *preLlr, double *appLlr, double **alpha, double **beta, double ***gamma, int sch);
double	LogMaxFunction(double x, double y);
void	TurboProcessor(LDPC &ldpc, int ***trellis, double *tempLlr, double **refLlr, double *rxLlr, double *preLlr, double **appLlr, double **alpha, double **beta, double ***gamma);

void	UpSampling(double ****tx, double ****preTx);
double	SquareRootRaisedCosine(double m);
void	SRRCGeneration(double *srrc, double drift);
void	PulseShaping(double ****preTx, double ****tx, double *txFilter , double *drift);
void	ReceivingFilter(double ***rx, double **tempRx, double *rxFilter);

void	EnergyProfile(double ****h, double ****H);
void	MultipleAccessChannel(double stdDev, double ****h, double ****tx, double ***rx, double *drift, double ****H);

void	MLDT(LDPC &ldpc, double variance, double ****H, double ***postRx, double **app, double **appLlr, double **refLlr, double ****estimate);

void	FFT(double *pr, double *pi, int n, int k, double *fr, double *fi, int l, int il);

void	Clustering(double ***postRx, double **centroid, int **group, int *groupSize, double *distList, double *variation, double **softAssign, double variance, double ****estimate);
void	InitialSeeding(double ***postRx, double **centroid, int **group, int *groupSize, double *distList, int sch, int offset);
void	Grouping(double ***postRx, double **centroid, int **group, int *groupSize, int sch, int offset);
void	CentroidRenewal(double ***postRx, double **centroid, int **group, int *groupSize, double *variation, int sch);
void	CoefEstimation(double **centroid, double ****estimate, int sch, int offset);
void	UserIdentification(double ****H, double ****estimate);
bool	ConditionCheck(double **centroid, int *groupSize, int &count);
void	EMClustering(double ***postRx, double **centroid, double **softAssign, double variance, int sch, int offset);
void    printchannel(double ****h);
template <class T>
void Swaping(T **x, T **y)
{
	T *temp = *x;
	*x = *y;
	*y = temp;
}