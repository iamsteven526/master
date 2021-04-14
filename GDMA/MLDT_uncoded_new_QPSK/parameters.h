#pragma once
#include <string>
#include <vector>
//---------- system ----------

#define		NUM_USER				4									// number of users
#define		NUM_LEVEL				( 1 << NUM_USER )					// number of levels of superimposed signal

#define		DIFF_ENC				0									// differential encoding; 1: enable, 0: disable
#define		DIFF_RX_SCHEME			0									// 0: differential decoding, 1: BCJR algo.

#define		Scrambler				0

//---------- channel estimation ----------

#define		CE_SCHEME				0									// 0: cluster-based, 1: ideal

#define		INI_METHOD				2									// 0: LBG, 1: k-means++, 2: modified k-means++

#define		DELTA					0.001								// spliting distance for LBG algo.

#define		GROUP_SIZE				NUM_LEVEL 						    // number of groups to be classified

#define		EM_GMM					1									// EM clustering with Gaussian mixture model; 1: enable, 0: disable

#define		PROPOSAL				2									// 0: original,1: proposal-1, 2: proposal-2

//---------- synchronization ----------

#define		SYNCHRONOUS				1									// 1: synchronous, 0: unknown timing drifts

#define		UP_RATE					10									// up-sampling rate

#define		TRUNCATION				6									// truncation of shaping filter
#define		ROLL_OFF				0.22								// roll-off factor

#define		DRIFT_RANGE				0.6									// drifting range

#define		MODIFIED_LLR			0									// 1: enable, 0: disable
#define		MODIFIER				1e-6								
#define		SPREAD_LEN	     	    7
#define		Power_Ratio		    	0.01
#define     Ndec					10
#define		SAMPLE_NUM				128

//---3 state [2,1]
#if SPREAD_LEN == 7

#endif

//---------- collision state-----

#define COLLISION					1									// 1:  Simultaneouly Transmit, 0: Transmit with known time drift 


//---------- simulator ----------

#define		NUMERIC_LIMIT			1e-100
#define		LLR_LIMIT				230

#define		BLOCK_NUM				600000						// number of blocks to be simulated
#define		BLOCK_LEN				1024							// number of symbols in a block

#define		SNR_NUM					9									// number of SNR points to be simulated
#define		SNR_START				30									// in dB
#define		SNR_STEP			    5									// in dB

#define		HARD(x)					( (x) > 0 ? 0 : 1 )

//--------------------

void	Transmitter(int **data, double *preTx, double **tx, double *txFilter, double *rxFilter, double** pilot, double* prePilot, int *known_drift);
void	Modulator(int *data, double *tx, int known_drift);
void	Detector(int **data, double **appLlr, int ***trellis, double **alpha, double **beta, double ***gamma, long double &error, int *known_drift);
void	DiffDecoding(double *appLlr, int known_drift);

void	UpSampling(double *tx, double *preTx, double* polit, double* prePilot);
double	SquareRootRaisedCosine(double m);
void	SRRCGeneration(double *srrc, double drift);
void	PulseShaping(double *preTx, double *tx, double *txFilter, double *rxFilter, double* prePilot, double* pilot);
void	ReceivingFilter(double **rx, double **postRx, double *rxFilter, double *txFilter, double** pilot, double **chCoef);

void	TrellisConstruction(int ***trellis);
void	BCJR(int ***trellis, double *rxLlr, double *appLlr, double **alpha, double **beta, double ***gamma);
double	LogMaxFunction(double x, double y);

void	EnergyProfile(double **chCoef);
void	MultipleAccessChannel(double stdDev, double **chCoef, double **tx, double **rx, double **pilot, int *known_drift);

void	MLDT(double variance, double **chCoef, double **rx, double **app, double **appLlr, double **estimate, int *known_drift);

void	Clustering(double **rx, double **centroid, int **group, int *groupSize, double *distList, double *variation, double **softAssign, double variance, double **estimate, long double &itCount, double **chCeof, int *known_drift);
void	InitialSeeding(double **rx, double **centroid, int **group, int *groupSize, double *distList, double variance, int CLUSTER_USER, int CLUSTER_GROUP, int CLUSTER_SIZE);
void	Grouping(double **rx, double **centroid, int **group, int *groupSize, int CLUSTER_USER, int CLUSTER_GROUP, int CLUSTER_SIZE);
double	EuclideanDistance(double *x, double *y);
void	CentroidRenewal(double **rx, double **centroid, int **group, int *groupSize, double *variation, int CLUSTER_GROUP);
void	CentroidRenewal2(double** rx, double** centroid, int** group, int* groupSize);
void	CoefEstimation(double **centroid, double **estimate, double variance, bool & reset);
void	MSEComparison(double **chCoef, double **estimate, double &mse, double **rx, double **centroid,int*known_drift);
bool	ConditionCheck(double **centroid, int *groupSize, int &count, int CLUSTER_GROUP);
bool	ConditionCheck1(double** centroid, int* groupSize, double variance, double** rx, double** estimate, double** chCoef, int &count, int **group);
void	EMClustering(double **rx, double **centroid, double **softAssign, double variance, int CLUSTER_LEN, int CLUSTER_GROUP);
void	printdata(double** rx, double** es, double** ideal, double** centroid);
void	printinitial(double** rx, double** centroid,double **chCoef,int cluster_size);
void	CentroidMSEComparison(double** chCoef, double** estimate, double& centroid_mse, double** rx, double** centroid);
void	Check_initial(double** rx, double** centroid, int** group, int* groupSize, double variance);
bool	Check_initial_modified(double** rx, double** centroid, int** group, int* groupSize, double variance, double **estimate, double **chCoef);
void	print_shape(double* txfilter, std::string name);
void	Inner_K_means(int min, int max,int sec_min, int* groupSize, int** group,double **centroid,double **rx);
std::vector<std::vector<double>> pair_seq(std::vector<std::vector<double>> centroid, int num_level, double& mse, std::vector<double> group_centroid, std::vector<int> pair_num);
std::vector<std::vector<double>> CoefEstimation_(double** centroid, double** estimate, double variance, int CLUSTER_GROUP, int GROUP_USER);
void	printdata2(double** rx, std::vector<std::vector<double>> es, double** chCoef, double** centroid, int CLUSTER_LEN, int CLUSTER_GROUP, int CLUSTER_USER);
void	scrambler(int* data);
std::vector<int> descrambler(double* data);

template <class T>
void Swaping(T **x, T **y)
{
	T *temp = *x;
	*x = *y;
	*y = temp;
}