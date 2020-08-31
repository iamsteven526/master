#ifndef _SETTING_
#define _SETTING_

#include<stdio.h>
#include<stdint.h>
#include<iostream>
#include<fstream>
#include<cmath>
#include<time.h>
#include<vector>
#include<complex>
#include<algorithm>
#include<math.h>
#include <cstring>

using namespace std;
typedef complex<double>	complexNum;

#define RadioFrameNum				800
#define SymbolPerFrame				140

#define LAYER						10									// the number of subcarriers (1024)
#define N							int(pow(2,LAYER))					
#define MOD_BIT						4									// 16-QAM

#define CHANNEL_TYPE				2									// 0:AWGN, 1: slow fading, 2: fast fading
#define	ALPHA						double(0.8)							// exponentially decay parameter(slow fading)					

#define SNR_dB_START				30									// Eb/No
#define INCREMENT					5
#define SNR_dB_END					30

#define L							6									// tap number
#define G							(L-1)								// CP length
#define D							10									// ICI depth
#define Ls							8									// list size
#define S							64									// segment length
#define UE_SPEED					540									// km/h 
#define CARRIER_FREQ				(3*pow(10,9))						// GHz 

#define R							1
#define DEVISION					1
#define H_FILE						"H_4096_2048_z128_0635.txt"
#define TANH_BOUND					50
#define ITERATION					50

#define DELRA_FREQ					(15*1000)							// subcarrier spacing (kHz)
#define Fs							(N*DELRA_FREQ)						// chip sampling rate (kHz)
#define CHIP_TIME					(1./Fs)
#define SYMBOL_DURATION				((N+G)*CHIP_TIME)

#define Es							float(1)
#define MIN_BLK_ERR					100
#define MIN_BIT_ERR					1000

template<class T>
T** CreateMatrix(int row, int col, T par)
{
	T **p = new T *[row];
	for (int i = 0; i < row; i++)
	{
		p[i] = new T[col];
		memset(p[i], 0, col * sizeof(T));
	}
	return p;
}

template<class T>
void DeleteMatrix(T **matrix)
{
	int row = _msize(matrix) / sizeof(matrix[0]);
	for (int i = 0; i < row; i++)
		delete[] matrix[i];
	delete[] matrix;
}

#endif // !_SETTING_

