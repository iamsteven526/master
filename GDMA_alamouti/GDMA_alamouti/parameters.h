#pragma once

//---------- system ----------

#define		NUM_USER				2													// number of users

#define		NUM_TX					2													// number of transmit antennas
#define		NUM_LEVEL				( 1 << ( NUM_USER * NUM_TX - 1))					// number of levels of superimposed signal

//---------- simulator ----------

#define		NUMERIC_LIMIT			1e-100
#define		LLR_LIMIT				230

#define		BLOCK_NUM				10000000											// number of blocks to be simulated

#define		SNR_NUM					9													// number of SNR points to be simulated
#define		SNR_START				0 													// in dB
#define		SNR_STEP				5													// in dB

#define		HARD(x)					( (x) > 0 ? 0 : 1 )

//--------------------

void	AlamoutiEncoder(int **data, double ***tx);
void	SignalCombiner(double ***chCoef, double **rx, double ***postRx);
void	ComplexConjugate(double *x, double *y);
void	ComplexMultiplication(double *x1, double *x2, double *y);
void	SuperLevelSpecification(double ***chCoef, double ****supLevel);
void	Detector(int **data, double **appLlr, long double &error);

void	EnergyProfile(double ***chCoef);
void	MultipleAccessChannel(double stdDev, double ***chCoef, double ***tx, double **rx);

void	MLDT(double variance, double ***chCoef, double ****supLevel, double ***postRx, double ***app, double **appLlr);