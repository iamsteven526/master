#pragma once
#include<vector>
//---------- system ----------

#define		NUM_USER				2									// number of users
#define		Qm						2									// QPSK
#define		NUM_LEVEL				( 1 << ( NUM_USER * Qm ) )			// number of levels of superimposed signal

#define		DIFF_ENC				0								// differential encoding; 1: enable, 0: disable

//---------- channel estimation setting ----------

#define		CE_METHOD				1									// 0: cluster-based; 1: ideal

#define		INI_METHOD				1									//  1: k-means++; 2: k-means++ modify; 

#define		DELTA					0.001								// spliting distance for LBG algo.
//#define		DELTA					0.05
#define		GROUP_SIZE				NUM_LEVEL							// number of groups to be classified

#define		EM_GMM					0									// EM clustering with Gaussian mixture models, 1: enable, 0: disable

//---------- simulator ----------

#define		NUMERIC_LIMIT			1e-100
#define		LLR_LIMIT				230

#define		BLOCK_NUM				1000000							// number of blocks to be simulated
#define		BLOCK_LEN				100								// number of symbols in a block
#define		Cluster_Len				100

#define		SNR_NUM					9									// number of SNR points to be simulated
#define		SNR_START			    0									// in dB
#define		SNR_STEP				5									// in dB

#define		HARD(x)					( (x) > 0 ? 0 : 1 )

//--------------------

void	Transmitter(int **data, double **tx);
void	Modulator(int *data, double *tx);
void	Detecter(int **data, double **appLlr, long double &error);
void	DiffDecoding(double *appLlr);

void	EnergyProfile(double **chCoef);
void	MultipleAccessChannel(double stdDev, double **chCoef, double **tx, double **rx);

void	MLDT(double variance, double **chCoef, double **rx, double *app, double **appLlr,double **estimate);

void	KmeansClustering(double** rx, double** centroid, int** group, int* groupSize, double* distList, double* variation, double** softAssign, double variance, double** estimate, long double& itCount, double** chCoef);
void	InitialSeeding(double** rx, double** centroid, int** group, int* groupSize, double* distList);
void	Grouping(double** rx, double** centroid, int** group, int* groupSize);
double	EuclideanDistance(double* x, double* y);
void	CentroidRenewal(double** rx, double** centroid, int** group, int* groupSize, double* variation);
void	CoefEstimation(double** centroid, double** estimate, double** rx, double** chCoef);
void	MSEComparison(double** chCoef, double** estimate, double& mse, double** rx, double** centroid);
bool	ConditionCheck(double** centroid, int* groupSize, int& count);
void	EMClustering(double** rx, double** centroid, double** softAssign, double variance, double** chCoef, double** estimate);
void	printdata(double** rx, double** es, double** ideal, double** centroid);
void	printinitial(double** rx, double** centroid);
void	CentroidMSEComparison(double** chCoef, double** estimate, double& centroid_mse, double** rx, double** centroid);
void    QPSKcheck(double ** centroid,int **pair,std::vector<std::vector<double>> &reg,double **estimate);
void	Check_initial(double** rx, double** centroid, int** group, int* groupSize);
void	Generate_square_sequence(std::vector<std::vector<double>> &centroid_per, int** pair, double** centroid, int pair_size,int sequence_num);
bool	Check_initial_modified(double** rx, double** centroid, int** group, int* groupSize, double variance, double** estimate, double** chCoef);
void	Inner_K_means(int min, int max, int sec_min, int* groupSize, int** group, double** centroid, double** rx);
bool    ConditionCheck1(double** centroid, int* groupSize, double variance, double** rx, double** estimate, double** chCoef, int& count, int** group);
void	CentroidRenewal2(double** rx, double** centroid, int** group, int* groupSize);
template <class T>
void Swaping(T** x, T** y)
{
	T* temp = *x;
	*x = *y;
	*y = temp;
}