#pragma once
#include<vector>
#include<string>
//---------- system ----------

#define		NUM_USER				3									// number of users
#define		Qm						2									// QPSK
#define		NUM_LEVEL				( 1 << ( NUM_USER * Qm ) )			// number of levels of superimposed signal

#define		DIFF_ENC				0								// differential encoding; 1: enable, 0: disable
#define		DIFF_RX_SCHEME			0									// 0: differential decoding, 1: BCJR algo.

#define		Scrambler				0
//-----------scma mode------//now only for perfectCSI

#define		SCMA_SOURCE				4									//default:1
#define		SCMA_USER_SOURCE		2									//sources for a user, default:1 




//---------- channel estimation setting ----------

#define		CE_METHOD				1									// 0: cluster-based; 1: ideal

#define		INI_METHOD				2									//  1: k-means++; 2: k-means++ modify; 

#define		DELTA					0.001								// spliting distance for LBG algo.
//#define		DELTA					0.05
#define		GROUP_SIZE				NUM_LEVEL							// number of groups to be classified

#define		EM_GMM					1									// EM clustering with Gaussian mixture models, 1: enable, 0: disable

#define		PROPOSAL				1									// 0: original,1: proposal-1, 2: proposal-2

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

//---------- collision state-----

#define COLLISION					1									// 1:  Simultaneouly Transmit, 0: Transmit with known time drift 

//---------- simulator ----------

#define		NUMERIC_LIMIT			1e-100
#define		LLR_LIMIT				230

#define		BLOCK_NUM				14000				// number of blocks to be simulated
#define		BLOCK_LEN				100							// number of symbols in a block
#define		Cluster_Len				100

#define		SNR_NUM					3									// number of SNR points to be simulated
#define		SNR_START			    10							// in dB
#define		SNR_STEP				5									// in dB

#define		HARD(x)					( (x) > 0 ? 0 : 1 )

//--------------------



void	Transmitter(int **data, double **tx);
void	Modulator(int *data, double *tx);
void	Detecter(int **data, double **appLlr, long double &error);
void	DiffDecoding(double *appLlr);

void	EnergyProfile(double ***chCoef);
void	MultipleAccessChannel(double stdDev, double ***chCoef, double*** chCoef2, double **tx, double **rx, int index, int scma_matrix[SCMA_SOURCE][SCMA_SOURCE * NUM_USER / SCMA_USER_SOURCE], int coef_idx[SCMA_SOURCE][SCMA_SOURCE * NUM_USER / SCMA_USER_SOURCE]);

void	fourwayMLDT(double variance, double*** chCoef, double*** chCoef2, double*** rx, double** app, double** appLlr, double*** finalestimate,double*** finalestimate2, int scma_matrix[SCMA_SOURCE][SCMA_SOURCE * NUM_USER / SCMA_USER_SOURCE], int coef_idx[SCMA_SOURCE][SCMA_SOURCE * NUM_USER / SCMA_USER_SOURCE]);
void	MLDT(double variance, double **chCoef, double** chCoef2, double **rx, double *app, double** appLlr, int i);
void	MLDTLLR(int i, double** app, double** appLlr);
void	normalization(double* app);


/*
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
*/











/*
void	Transmitter(int **data, double *preTx, double **tx, double *txFilter, double *rxFilter, double** pilot, double* prePilot, int *known_drift);
//void	Transmitter(int** data, double* preTx, double** tx, double* txFilter,double **pilot);
void	Modulator(int *data, double *tx, int known_drift);
void	Detector(int **data, double **appLlr, int ***trellis, double **alpha, double **beta, double ***gamma, long double &error, int *known_drift);
void	DiffDecoding(double *appLlr, int known_drift);
*/

void	UpSampling(double *tx, double *preTx, double* polit, double* prePilot);
double	SquareRootRaisedCosine(double m);
void	SRRCGeneration(double *srrc, double drift);
void	PulseShaping(double *preTx, double *tx, double *txFilter, double *rxFilter, double* prePilot, double* pilot);
void	ReceivingFilter(double **rx, double **postRx, double *rxFilter, double *txFilter, double** pilot, double **chCoef);


void	TrellisConstruction(int ***trellis);
void	BCJR(int ***trellis, double *rxLlr, double *appLlr, double **alpha, double **beta, double ***gamma);
double	LogMaxFunction(double x, double y);


/*
void	EnergyProfile(double **chCoef);
void	MultipleAccessChannel(double stdDev, double **chCoef, double **tx, double **rx, double **pilot, int *known_drift);

void	MLDT(double variance, double **chCoef, double **rx, double **app, double **appLlr, double **estimate, int *known_drift);
*/


void	Clustering(double **rx, double **centroid, int **group, int *groupSize, double *distList, double *variation, double **softAssign, double variance, double **estimate, long double &itCount, double ***chCeof, int *known_drift);
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
void    SCMA_MSEComparison(double*** chCoef,double*** chCoef2,double** estimate,double*** finalestimate,double*** finalestimate2,int a[SCMA_SOURCE * NUM_USER / SCMA_USER_SOURCE], int b[SCMA_SOURCE * NUM_USER / SCMA_USER_SOURCE]);
/*
void	scrambler(int* data);
std::vector<int> descrambler(double* data);
*/


void	TrellisConstruction(int ***trellis);
void	BCJR(int ***trellis, double *rxLlr, double *appLlr, double **alpha, double **beta, double ***gamma);
double	LogMaxFunction(double x, double y);



template <class T>
void Swaping(T** x, T** y)
{
	T* temp = *x;
	*x = *y;
	*y = temp;
}