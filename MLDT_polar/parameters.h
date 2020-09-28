#pragma once

#include <vector>
#include "polar.h"

//---------- Polar code ----------

#define		BCT_LAYER				7								    // number of layers in channel polarization
#define		DATA_LEN				64								    // data length
#define		CODE_LEN				( 1 << BCT_LAYER )					// codeword length

#define		CRC_LEN					0									// number of bits in CRC
#define		LIST_SIZE				32									// list size 

#define		PUNCTURING				0									// 1: enable, 0: disable
#define		PUNCTURED_BIT			16									// number of bits to be punctured

#define		PCC_METHOD				0									// 0: Bhattacharyya bound
																		// 1: Capacity bound
																		// 2: Piecewise linear approximation
																		// 3: Gaussain approximation
																		// 4: 3GPP; fixed lengh of 1024
																		// 5: NBC; fixed lengh of 1024	

#define		DESIGHED_SNR			1									// SNR per transmitted symbol
#define		DECODING_TYPE			1									// 1: SCL,  0: BP
#define		JOINT					0								    // 1: JCD,  0: SCD
#define		Iteration				50
//---------- system ----------
#define		Racian					0												// 1 : Rician
																					// 0 : Rayleight
#define		K						3
#define     phase					11

#define		NUM_USER				3									// number of users
#define		NUM_LEVEL				( 1 << NUM_USER )					// number of levels of superimposed signal

#define		NBC						1									// noncoherent block coding; 1: enable, 0: disable 

//---------- channel estimation ----------

#define		CE_SCHEME				0									// 0: cluster-based, 1: ideal

#define		INI_METHOD				2									// 0: LBG, 1: k-means++, 2: modified k-means++

#define		GROUP_SIZE				NUM_LEVEL							// number of groups to be classified

#define		EM_GMM					1									// EM clustering with Gaussian mixture model; 1: enable, 0: disable

#define		PROPOSAL				1									// 0: original,1: proposal-1, 2: proposal-2

#define		DELTA					0.001								// spliting distance for LBG algo.
//---------- synchronization ----------

#define		SYNCHRONOUS				1									// 1: synchronous, 0: unknown timing drifts

#define		UP_RATE					4									// up-sampling rate

#define		TRUNCATION				6									// truncation of shaping filter	 							
#define		ROLL_OFF				0.22								// roll-off factor

#define		DRIFT_RANGE				0.1									// drifting range

#define		MODIFIED_LLR			0									// 1: enable, 0: disable
#define		MODIFIER				1e-6

#define		COLLISION				1									// 1. transmit simultanously , 2. known time drift
//---------- simulator ----------

#define		NUMERIC_LIMIT			1e-100
#define		LLR_LIMIT				230

#define		BLOCK_NUM				1000000								// number of blocks to be simulated

#define		SNR_NUM					7									// number of SNR points to be simulated
#define		SNR_START				15									// in dB
#define		SNR_STEP				5									// in dB

#define		HARD(x)					( (x) > 0 ? 0 : 1 )

// ----------channel model----------
#define		FADING_TYPE				0												// 0: block fading 1: doppler fading 2: Qusic-static
#define		FFT_POINT				16

#if FADING_TYPE == 0
	#define		FADING_SIZE				CODE_LEN
#elif FADING_TYPE ==1
	#define		FADING_SIZE				1
#elif FADING_TYPE == 2
	#define		FADING_SIZE			    64
#endif

#define		CARRIER_FREQ			2												// carrier frequency (GHz)
#define		CHIP_RATE				1.92 //30.72									// sampling frequency (MHz)

#define		UE_SPEED				300												// mobile speed (km/hour)


//--------------------

void	Transmitter(PolarCode &polar, std::vector<std::vector<uint8_t>> &data, std::vector<std::vector<uint8_t>> &codeword, double *preTx, double **tx, double *txFilter);
void	Modulator(std::vector<uint8_t> &codeword, double *tx);
void	Detector(PolarCode &polar, std::vector<std::vector<uint8_t>> &data, std::vector<std::vector<double>> &appLlr, long double &errBlock, long double &errBit, double **app, double ***ch, double *llr_spa);

void	UpSampling(double *tx, double *preTx);
double	SquareRootRaisedCosine(double m);
void	SRRCGeneration(double *srrc, double drift);
void	PulseShaping(double *preTx, double *tx, double *txFilter);
void	ReceivingFilter(double **rx, double **postRx, double *rxFilter);

void	EnergyProfile(double ***chCoef);
void	MultipleAccessChannel(double stdDev, double ***chCoef, double **tx, double **rx);

void	MLDT(double variance, double ***chCoef, double **rx, double **app, std::vector<std::vector<double>> &appLlr, double **estimate);

void	Clustering(double** rx, double** centroid, int** group, int* groupSize, double* distList, double* variation, double** softAssign, double variance, double** estimate, double*** chCoef);
void	InitialSeeding(double** rx, double** centroid, int** group, int* groupSize, double* distList, double variance, int CLUSTER_USER, int CLUSTER_GROUP, int CLUSTER_SIZE);
void	Grouping(double** rx, double** centroid, int** group, int* groupSize, int CLUSTER_USER, int CLUSTER_GROUP, int CLUSTER_SIZE);
double	EuclideanDistance(double* x, double* y);
void	CentroidRenewal(double** rx, double** centroid, int** group, int* groupSize, double* variation, int CLUSTER_GROUP);
void	CentroidRenewal2(double** rx, double** centroid, int** group, int* groupSize);
void	CoefEstimation(double** centroid, double** estimate, double variance, bool& reset);
void	MSEComparison(double ***chCoef, double **estimate, double &mse);
bool	ConditionCheck(double** centroid, int* groupSize, int& count, int CLUSTER_GROUP);
void	EMClustering(double** rx, double** centroid, double** softAssign, double variance, int CLUSTER_LEN, int CLUSTER_GROUP);
void	Inner_K_means(int min, int max, int sec_min, int* groupSize, int** group, double** centroid, double** rx);
void	printinitial(double** rx, double** centroid, double*** chCoef);
std::vector<std::vector<double>> pair_seq(std::vector<std::vector<double>> centroid, int num_level, double &mse, std::vector<double> group_centroid, std::vector<int> pair_num);

template <class T>
void Swaping(T **x, T **y)
{
	T *temp = *x;
	*x = *y;
	*y = temp;
}