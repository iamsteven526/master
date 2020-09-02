#define _USE_MATH_DEFINES
#include <iostream>
#include "parameters.h"
using namespace std;

//20200211 add NUM_USER == 3 with CE_METHOD = 1
void MLDT(double variance, double **chCoef, double **rx, double *app, double **appLlr,double **estimate) //chCoef[userid][modulationlevel]
{
	double chCoefreg = 0, test = 0;
	double probzero = 0 ,probone = 0;
	if (CE_METHOD)
	{
		for (int i = 0; i < BLOCK_LEN; i++)
		{
			//---------- APPs ----------//real
			if (NUM_USER == 1)
			{
				app[0] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1])), 2.) / (2. * variance));
				app[1] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1])), 2.) / (2. * variance));
				app[2] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1])), 2.) / (2. * variance));
				app[3] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1])), 2.) / (2. * variance));
			}
			else if (NUM_USER == 2)
			{
				app[0] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[1] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[2] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[3] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[4] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[5] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[6] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[7] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[8] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[9] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[10] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[11] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[12] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[13] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[14] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[15] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
			}
			else if (NUM_USER >= 3)//20200211 bug....QQ
			{
				for (int appcounts = 0; appcounts < pow(2,(Qm * NUM_USER)); ++appcounts) 
				{
					for (int userid = 0; userid < NUM_USER; ++userid)
					{
						for (int qmlevel = 0; qmlevel < Qm; ++qmlevel) 
						{
							test = (Qm * NUM_USER) - userid * 2 - qmlevel;//base for real ex:( - + - + - +)(6 5 4 3 2 1)
							chCoefreg += pow(-1, test + 1) * pow(-1, appcounts / int(pow(2, test - 1))) * chCoef[userid][qmlevel];
						}	
					}
					//cout << chCoefreg << endl;//debug
					app[appcounts] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoefreg)), 2.) / (2. * variance));
					chCoefreg = 0;
				}
			}
			else
			{
				printf("\nPARAMETER SETTING IS WRONG\n");
				system("pause");
			}
			//---------- normalization ----------
			double temp = 0;
			for (int j = 0; j < NUM_LEVEL; j++)
			{
				if (app[j] < NUMERIC_LIMIT) app[j] = NUMERIC_LIMIT;
				temp += app[j];
			}
			for (int j = 0; j < NUM_LEVEL; j++)
			{
				app[j] /= temp;
				if (app[j] < NUMERIC_LIMIT) app[j] = NUMERIC_LIMIT;
			}
			//---------- APPs ----------//imag
			if (NUM_USER == 1)
			{
				app[0] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1])), 2.) / (2. * variance));
				app[1] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1])), 2.) / (2. * variance));
				app[2] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1])), 2.) / (2. * variance));
				app[3] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1])), 2.) / (2. * variance));
			}
			else if (NUM_USER == 2)
			{
				app[0] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[1] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[2] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[3] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[4] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[5] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[6] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[7] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[8] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[9] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[10] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[11] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[12] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[13] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[14] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[15] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
			}
			else if (NUM_USER >= 3)//20200211
			{
				for (int appcounts = 0; appcounts < pow(2, (Qm * NUM_USER)); ++appcounts)
				{
					for (int userid = 0; userid < NUM_USER; ++userid)
					{
						for (int qmlevel = 0; qmlevel < Qm; ++qmlevel)
						{
							test = (Qm * NUM_USER) - userid * 2 - (1 - qmlevel);//base for imag ex:( - - - - - -)(5 6 3 4 1 2)
							chCoefreg += (-1) * pow(-1, appcounts / int(pow(2, test - 1))) * chCoef[userid][qmlevel];
						}
					}
					app[appcounts] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoefreg)), 2.) / (2. * variance));
					chCoefreg = 0;
				}
			}

			//---------- normalization ----------
			temp = 0;
			for (int j = 0; j < NUM_LEVEL; j++)
			{
				if (app[j] < NUMERIC_LIMIT) app[j] = NUMERIC_LIMIT;
				temp += app[j];
			}
			for (int j = 0; j < NUM_LEVEL; j++)
			{
				app[j] /= temp;
				if (app[j] < NUMERIC_LIMIT) app[j] = NUMERIC_LIMIT;
			}
			//---------- a-posteriori LLRs ----------
			if (NUM_USER == 1)
			{
				if ((app[0] + app[1]) <= NUMERIC_LIMIT) appLlr[0][2 * i] = -LLR_LIMIT;
				else if ((app[2] + app[3]) <= NUMERIC_LIMIT) appLlr[0][2 * i] = LLR_LIMIT;
				else appLlr[0][2 * i] = log((app[0] + app[1]) / (app[2] + app[3]));
				if ((app[0] + app[2]) <= NUMERIC_LIMIT) appLlr[0][2 * i + 1] = -LLR_LIMIT;
				else if ((app[1] + app[3]) <= NUMERIC_LIMIT) appLlr[0][2 * i + 1] = LLR_LIMIT;
				else appLlr[0][2 * i + 1] = log((app[0] + app[2]) / (app[1] + app[3]));
			}
			else if (NUM_USER == 2)
			{
				//---------- user-A ----------
				if ((app[0] + app[1] + app[2] + app[3] + app[4] + app[5] + app[6] + app[7]) <= NUMERIC_LIMIT) appLlr[0][2 * i] = -LLR_LIMIT;
				else if ((app[8] + app[9] + app[10] + app[11] + app[12] + app[13] + app[14] + app[15]) <= NUMERIC_LIMIT) appLlr[0][2 * i] = LLR_LIMIT;
				else appLlr[0][2 * i] = log((app[0] + app[1] + app[2] + app[3] + app[4] + app[5] + app[6] + app[7]) / (app[8] + app[9] + app[10] + app[11] + app[12] + app[13] + app[14] + app[15]));
				if ((app[0] + app[1] + app[2] + app[3] + app[8] + app[9] + app[10] + app[11]) <= NUMERIC_LIMIT) appLlr[0][2 * i + 1] = -LLR_LIMIT;
				else if ((app[4] + app[5] + app[6] + app[7] + app[12] + app[13] + app[14] + app[15]) <= NUMERIC_LIMIT) appLlr[0][2 * i + 1] = LLR_LIMIT;
				else appLlr[0][2 * i + 1] = log((app[0] + app[1] + app[2] + app[3] + app[8] + app[9] + app[10] + app[11]) / (app[4] + app[5] + app[6] + app[7] + app[12] + app[13] + app[14] + app[15]));
				//---------- user-B ----------
				if ((app[0] + app[1] + app[4] + app[5] + app[8] + app[9] + app[12] + app[13]) <= NUMERIC_LIMIT) appLlr[1][2 * i] = -LLR_LIMIT;
				else if ((app[2] + app[3] + app[6] + app[7] + app[10] + app[11] + app[14] + app[15]) <= NUMERIC_LIMIT) appLlr[1][2 * i] = LLR_LIMIT;
				else appLlr[1][2 * i] = log((app[0] + app[1] + app[4] + app[5] + app[8] + app[9] + app[12] + app[13]) / (app[2] + app[3] + app[6] + app[7] + app[10] + app[11] + app[14] + app[15]));
				if ((app[0] + app[2] + app[4] + app[6] + app[8] + app[10] + app[12] + app[14]) <= NUMERIC_LIMIT) appLlr[1][2 * i + 1] = -LLR_LIMIT;
				else if ((app[1] + app[3] + app[5] + app[7] + app[9] + app[11] + app[13] + app[15]) <= NUMERIC_LIMIT) appLlr[1][2 * i + 1] = LLR_LIMIT;
				else appLlr[1][2 * i + 1] = log((app[0] + app[2] + app[4] + app[6] + app[8] + app[10] + app[12] + app[14]) / (app[1] + app[3] + app[5] + app[7] + app[9] + app[11] + app[13] + app[15]));
			}
			else if (NUM_USER >= 3)
			{
				for (int userid = 0; userid < NUM_USER; ++userid) {
					for (int bitnum = 0; bitnum < 2; ++bitnum) {
						for (int p = 0; p < pow(2, Qm * NUM_USER); ++p) {
							test = (Qm * NUM_USER) - userid * 2 - bitnum;//ex(6 5 4 3 2 1)
							if (pow(-1, p / int(pow(2, test - 1))) == 1) probzero += app[p];
							else probone += app[p];
						}
						if ( probzero <= NUMERIC_LIMIT) appLlr[userid][2 * i + bitnum] = -LLR_LIMIT;
						else if ( probone <= NUMERIC_LIMIT) appLlr[userid][2 * i + bitnum] = LLR_LIMIT;
						else appLlr[userid][2 * i + bitnum] = log(probzero / probone);
						//cout << probzero << endl;
						probzero = 0;
						probone = 0;
					}
				}
				//---------- user-userid ----------
			}
		}
	}
	else
	{
		for (int i = 0; i < NUM_USER; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				chCoef[i][j]=estimate[i][j];
			}
		}
		for (int i = 0; i < BLOCK_LEN; i++)
		{
			//---------- APPs ----------
			if (NUM_USER == 1)
			{
				app[0] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1])), 2.) / (2. * variance));
				app[1] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1])), 2.) / (2. * variance));
				app[2] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1])), 2.) / (2. * variance));
				app[3] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1])), 2.) / (2. * variance));
			}
			else if (NUM_USER == 2)
			{
				app[0] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[1] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[2] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[3] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[4] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[5] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[6] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[7] = exp(-pow((rx[i][0] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[8] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[9] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[10] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[11] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[12] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[13] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[14] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[15] = exp(-pow((rx[i][0] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
			}
			else
			{
				printf("\nPARAMETER SETTING IS WRONG\n");
				system("pause");
			}
			//---------- normalization ----------
			double temp = 0;
			for (int j = 0; j < NUM_LEVEL; j++)
			{
				if (app[j] < NUMERIC_LIMIT) app[j] = NUMERIC_LIMIT;
				temp += app[j];
			}
			for (int j = 0; j < NUM_LEVEL; j++)
			{
				app[j] /= temp;
				if (app[j] < NUMERIC_LIMIT) app[j] = NUMERIC_LIMIT;
			}
			//---------- APPs ----------
			if (NUM_USER == 1)
			{
				app[0] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1])), 2.) / (2. * variance));
				app[1] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1])), 2.) / (2. * variance));
				app[2] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1])), 2.) / (2. * variance));
				app[3] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1])), 2.) / (2. * variance));
			}
			else if (NUM_USER == 2)
			{
				app[0] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[1] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[2] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[3] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[4] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[5] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[6] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[7] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[8] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[9] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[10] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[11] *= exp(-pow((rx[i][1] - sqrt(1. / 2.) * (chCoef[0][0] - chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[12] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[13] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] - chCoef[1][1])), 2.) / (2. * variance));
				app[14] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] - chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
				app[15] *= exp(-pow((rx[i][1] + sqrt(1. / 2.) * (chCoef[0][0] + chCoef[0][1] + chCoef[1][0] + chCoef[1][1])), 2.) / (2. * variance));
			}
			//---------- normalization ----------
			temp = 0;
			for (int j = 0; j < NUM_LEVEL; j++)
			{
				if (app[j] < NUMERIC_LIMIT) app[j] = NUMERIC_LIMIT;
				temp += app[j];
			}
			for (int j = 0; j < NUM_LEVEL; j++)
			{
				app[j] /= temp;
				if (app[j] < NUMERIC_LIMIT) app[j] = NUMERIC_LIMIT;
			}
			//---------- a-posteriori LLRs ----------
			if (NUM_USER == 1)
			{
				if ((app[0] + app[1]) <= NUMERIC_LIMIT) appLlr[0][2 * i] = -LLR_LIMIT;
				else if ((app[2] + app[3]) <= NUMERIC_LIMIT) appLlr[0][2 * i] = LLR_LIMIT;
				else appLlr[0][2 * i] = log((app[0] + app[1]) / (app[2] + app[3]));
				if ((app[0] + app[2]) <= NUMERIC_LIMIT) appLlr[0][2 * i + 1] = -LLR_LIMIT;
				else if ((app[1] + app[3]) <= NUMERIC_LIMIT) appLlr[0][2 * i + 1] = LLR_LIMIT;
				else appLlr[0][2 * i + 1] = log((app[0] + app[2]) / (app[1] + app[3]));
			}
			else if (NUM_USER == 2)
			{
				//---------- user-A ----------
				if ((app[0] + app[1] + app[2] + app[3] + app[4] + app[5] + app[6] + app[7]) <= NUMERIC_LIMIT) appLlr[0][2 * i] = -LLR_LIMIT;
				else if ((app[8] + app[9] + app[10] + app[11] + app[12] + app[13] + app[14] + app[15]) <= NUMERIC_LIMIT) appLlr[0][2 * i] = LLR_LIMIT;
				else appLlr[0][2 * i] = log((app[0] + app[1] + app[2] + app[3] + app[4] + app[5] + app[6] + app[7]) / (app[8] + app[9] + app[10] + app[11] + app[12] + app[13] + app[14] + app[15]));
				if ((app[0] + app[1] + app[2] + app[3] + app[8] + app[9] + app[10] + app[11]) <= NUMERIC_LIMIT) appLlr[0][2 * i + 1] = -LLR_LIMIT;
				else if ((app[4] + app[5] + app[6] + app[7] + app[12] + app[13] + app[14] + app[15]) <= NUMERIC_LIMIT) appLlr[0][2 * i + 1] = LLR_LIMIT;
				else appLlr[0][2 * i + 1] = log((app[0] + app[1] + app[2] + app[3] + app[8] + app[9] + app[10] + app[11]) / (app[4] + app[5] + app[6] + app[7] + app[12] + app[13] + app[14] + app[15]));
				//---------- user-B ----------
				if ((app[0] + app[1] + app[4] + app[5] + app[8] + app[9] + app[12] + app[13]) <= NUMERIC_LIMIT) appLlr[1][2 * i] = -LLR_LIMIT;
				else if ((app[2] + app[3] + app[6] + app[7] + app[10] + app[11] + app[14] + app[15]) <= NUMERIC_LIMIT) appLlr[1][2 * i] = LLR_LIMIT;
				else appLlr[1][2 * i] = log((app[0] + app[1] + app[4] + app[5] + app[8] + app[9] + app[12] + app[13]) / (app[2] + app[3] + app[6] + app[7] + app[10] + app[11] + app[14] + app[15]));
				if ((app[0] + app[2] + app[4] + app[6] + app[8] + app[10] + app[12] + app[14]) <= NUMERIC_LIMIT) appLlr[1][2 * i + 1] = -LLR_LIMIT;
				else if ((app[1] + app[3] + app[5] + app[7] + app[9] + app[11] + app[13] + app[15]) <= NUMERIC_LIMIT) appLlr[1][2 * i + 1] = LLR_LIMIT;
				else appLlr[1][2 * i + 1] = log((app[0] + app[2] + app[4] + app[6] + app[8] + app[10] + app[12] + app[14]) / (app[1] + app[3] + app[5] + app[7] + app[9] + app[11] + app[13] + app[15]));
			}
		}
	}
}