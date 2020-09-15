#include "parameters.h"

#if DEGREE_SET == 311 // systematic, designed SNR = 1dB
#define STEP_SIZE 1
#define DEG_NUM 5
int step[STEP_SIZE] = { 0 };
int deg[DEG_NUM] = { 0, 1, 13, 14, 100 };
double deg_pdf[STEP_SIZE][DEG_NUM] = { 0, 0, 0.7128, 0.2783, 0.0089 };
#endif

#if DEGREE_SET == 312 // systematic, designed SNR = 3dB
#define STEP_SIZE 1
#define DEG_NUM 4
int step[STEP_SIZE] = { 0 };
int deg[DEG_NUM] = { 0, 1, 30, 31 };
double deg_pdf[STEP_SIZE][DEG_NUM] = { 0, 0, 0.5420, 0.4580 };
#endif

#if DEGREE_SET == 313 // systematic, designed SNR = -1dB
#define	STEP_SIZE 1
#define	DEG_NUM	5
int step[STEP_SIZE] = { 0 };
int deg[DEG_NUM] = { 0, 1, 7, 8, 100 };
double deg_pdf[STEP_SIZE][DEG_NUM] = { 0, 0, 0.1585, 0.8412, 0.0002 };
#endif

#if DEGREE_SET == 4717 // systematic; adaptive; Hui-Ming Wang, master thesis, 2018
#define STEP_SIZE 3
#define DEG_NUM 8
#if LDPC_H_COL == 2000
int step[] = { 2250, 2450 };
#elif LDPC_H_COL == 4000
int step[] = { 4490, 4890 };
#endif
int deg[DEG_NUM] = { 0, 1, 6, 7, 12, 13, 19, 20 };
double deg_pdf[STEP_SIZE][DEG_NUM] = { \
{0,		0,		0,		0,		0,		0, 0.0572, 0.9428},\
{0,		0,		0,		0, 0.8946, 0.1054,		0,		0},\
{0,		0, 0.1282, 0.8718,		0,		0,		0,		0} };
#endif

#if DEGREE_SET == 6221 // non-systematic; adaptive, s = 6
#define STEP_SIZE 7
#define DEG_NUM 18
#if LDPC_H_COL == 2000
int step[] = { 2235, 2850, 3800, 4560, 5570, 9500 };
#elif LDPC_H_COL == 10000
int step[] = { 11177, 14250, 19000, 22800, 28500, 47500 };
#endif
int deg[DEG_NUM] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 20, 45, 50, 55, 60, 65, 70 };
double deg_pdf[STEP_SIZE][DEG_NUM] = { \
{0, 0.0054, 0.4850, 0.1701, 0.1261,      0, 0.0238, 0.0940,      0,      0, 0.0630,      0, 0.0051, 0.0276,      0,      0,      0,      0},\
{0, 0.0110, 0.4099, 0.2063, 0.1020,      0, 0.0859, 0.0657,      0,      0, 0.0834,      0,      0,      0, 0.0334, 0.0023,      0,      0},\
{0, 0.0142, 0.3942, 0.2262, 0.0737,      0, 0.1224, 0.0447,      0,      0, 0.0794, 0.0068,      0, 0.0088, 0.0296, 0.0001,      0,      0},\
{0, 0.0277, 0.3329, 0.2526, 0.0648,      0, 0.1042, 0.0808,      0,      0, 0.0747, 0.0191,      0,      0, 0.0236, 0.0196,      0,      0},\
{0, 0.0442, 0.3202, 0.2595, 0.0482,      0, 0.1367, 0.0512,      0,      0, 0.0746, 0.0271,      0,      0, 0.0001, 0.0001, 0.0260, 0.0120},\
{0, 0.0509, 0.3519, 0.2369, 0.0447, 0.0932,      0, 0.0001, 0.1205, 0.0002, 0.0075, 0.0623,      0,      0,      0,      0, 0.0317,      0},\
{0, 0.0313, 0.3859, 0.2221, 0.0725, 0.0373, 0.0481, 0.0435, 0.0482, 0.0001, 0.0461, 0.0304, 0.0012, 0.0073, 0.0070, 0.0017, 0.0158, 0.0014} };
#endif

#if DEGREE_SET == 1 // non-systematic; MET; designed SNR = 0dB, jointly decoding
#define	STEP_SIZE 1
#define	DEG_NUM	11
int step[STEP_SIZE] = { 0 };
int deg[DEG_NUM] = { 0, 1, 2, 3, 4, 5, 7, 12, 14, 31, 50 };
double deg_pdf[STEP_SIZE][DEG_NUM] = { 0, 0.0188, 0.4031, 0.2758, 0.1005, 0.0406, 0.0425, 0.0517, 0.0621, 0.0039, 0.0010 };
#endif