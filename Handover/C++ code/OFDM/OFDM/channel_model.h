#ifndef CHANNEL_MODEL
#define CHANNEL_MODEL
#define _USE_MATH_DEFINES
#include"setting.h"
#include <random>
using namespace std;
random_device seed;
mt19937 generator(seed());
uniform_real_distribution<double> uniform(-0.5, 0.5);
normal_distribution<double> normal(0, 1);

class CHANNEL
{
public:

	vector<complexNum> ExponentiallyDecayChannel(int pathNum, double a)
	{
		vector<complexNum> complexGuassionImpulseResponse(pathNum);
		double ampitude = sqrt((1 - a) / (1 - pow(a, pathNum)));

		for (int i = 0; i < pathNum; i++)
			complexGuassionImpulseResponse.at(i) = ampitude * pow(a, float(i) / 2) * complexNum(sqrt(0.5)*normal(generator), sqrt(0.5)*normal(generator));

		return complexGuassionImpulseResponse;
	}

	vector<complexNum> QuasiStaticChannel(vector<complexNum> inputSequence, vector<complexNum> complexGuassionImpulseResponse, vector<double> &symbolPower, int symbolNum)
	{
		if (symbolNum % SymbolPerFrame == 0) {
			double power = 0;
			for (int l = 0; l < complexGuassionImpulseResponse.size(); l++)
				power += norm(complexGuassionImpulseResponse.at(l));
			symbolPower.push_back(power);
		}

		int initialLength = inputSequence.size(); 
		int convolutionLength = initialLength + complexGuassionImpulseResponse.size() - 1;  

		vector<complexNum> channelConvolationoutput(convolutionLength);
		complexNum inputPoint, channelGain;
		for (int i = 0; i < channelConvolationoutput.size(); i++)
		{
			for (int j = 0; j <= i; j++)
			{
				(j >= initialLength) ? inputPoint = 0 : inputPoint = inputSequence.at(j);
				((i - j) >= complexGuassionImpulseResponse.size()) ? channelGain = 0 : channelGain = complexGuassionImpulseResponse.at(i - j);
				channelConvolationoutput.at(i) += inputPoint * channelGain;
			}
		}

		return channelConvolationoutput;
	}

	vector<vector<complexNum>> TimeVaryingRayleighFadingChannel(vector<complexNum> fixedTimeCIR, int OFDMSymbolIndex)
	{
		vector<vector<complexNum>> timeVaryingChannelImpulseResponse(fixedTimeCIR.size(), vector<complexNum>((N + G + L - 1), 0));	// row = tap = fixedTimeCIR.size(), col = samples/symbol = N+L-1;
		
		double chip_period = 1. / Fs;
		double dopplerSpread = CARRIER_FREQ * (UE_SPEED*(1000. / 3600)) / (3. * pow(10, 8));
		
		int M = 16;		// Jakes' model parameter (resolution, fast fading)

		if (OFDMSymbolIndex == 1)
		{
			sita.resize(L);
			JakesAlpha.resize(0);

			vector<double> subAlpha(M);
			for (int tap = 0; tap < L; tap++)
			{
				sita.at(tap) = 2 * M_PI*uniform(generator);
				for (int m = 0; m < M; m++) 
					subAlpha.at(m) = (2* M_PI*(m + 1) - M_PI + sita.at(tap)) / (4 * M);
				
				JakesAlpha.push_back(subAlpha);
			}

			initialPhase1.resize(M);
			initialPhase2.resize(M);
			for (int m = 0; m < M; m++)
			{
				initialPhase1.at(m) = 2 * M_PI*uniform(generator);
				initialPhase2.at(m) = 2 * M_PI*uniform(generator);
			}
		}

		for (int tap = 0; tap < fixedTimeCIR.size(); tap++)
		{
			double startTime;
			if (OFDMSymbolIndex == 0)
				startTime = SYMBOL_DURATION * (SymbolPerFrame - 1);
			else
				startTime = SYMBOL_DURATION * (OFDMSymbolIndex - 1);

			for (int idx = 0; idx < timeVaryingChannelImpulseResponse[tap].size(); idx++)
			{
				double timeInstant = startTime + idx * chip_period;
				double Inpahse = 0, Quadrature = 0;

				for (int m = 0; m < M; m++)
				{
					Inpahse += cos(2 * M_PI*dopplerSpread*timeInstant*cos(JakesAlpha[tap].at(m)) + initialPhase1.at(m));
					Quadrature += cos(2 * M_PI*dopplerSpread*timeInstant*sin(JakesAlpha[tap].at(m)) + initialPhase2.at(m));
				}
				timeVaryingChannelImpulseResponse[tap].at(idx) = fixedTimeCIR.at(tap)*complexNum(sqrt(1. / M)*Inpahse, sqrt(1. / M)*Quadrature);

			}
		}

		return timeVaryingChannelImpulseResponse;
	}

	vector<complexNum> FastFadingChannel(vector<complexNum> inputSequence, vector<vector<complexNum>> TimeVaryingCIR, vector<double> &symbolPower)
	{
		int taps = TimeVaryingCIR.size();
		int timeidx = TimeVaryingCIR.at(0).size();
		double power = 0;

		for (int r = 0; r < taps; r++)
		{
			for (int c = 0; c < timeidx; c++)
				power += norm(TimeVaryingCIR.at(r).at(c));
		}
		symbolPower.push_back(power / timeidx);

		int initialLength = inputSequence.size();	

		vector<complexNum> channelConvolationoutput(timeidx, 0);

		for (int n = 0; n < timeidx; n++)
		{
			for (int l = 0; l < L; l++)
			{
				complexNum input = 0;
				((n - l) < 0 || (n - l) >= (N + G)) ? input = 0 : input = inputSequence.at(n - l);
				channelConvolationoutput.at(n) += TimeVaryingCIR.at(l).at(n)*input;
			}
		}

		return channelConvolationoutput;
	}

	vector<complexNum> AWGNNoise(int length, double	N0)
	{

		vector<complexNum> complexGuassionNoise(length);
		for (int i = 0; i < complexGuassionNoise.size(); i++)
			complexGuassionNoise.at(i) = complexNum(sqrt(N0 / 2)*normal(generator), sqrt(N0 / 2)*normal(generator));

		return complexGuassionNoise;
	}

private:
	vector<double> sita, initialPhase1, initialPhase2;

	vector<vector<double>> JakesAlpha;

};


#endif