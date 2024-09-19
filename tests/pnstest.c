// MIT License
// Copyright (c) 2024 Lauri Lorenzo Fiestas
// https://github.com/PrinssiFiestas/DSP/blob/main/LICENCE

#include <dspapi.h>
#include <chelpers.h>

#include <stdbool.h>

DSP_EXPORT void*   host=null;
DSP_EXPORT HostPrintFunc* hostPrint=null;
DSP_EXPORT double sampleRate = 0;
DSP_EXPORT uint audioInputsCount = 0;

//-------------------------------------------------------

DSP_EXPORT string const name="Kokeiluu";
DSP_EXPORT string description="Testailuu";

enum {
    PARAM_FREQ,
    PARAM_Q,
    PARAM_VOL,
    PARAM_GAIN,
    PARAMS_LENGTH
};

double params[PARAMS_LENGTH];
DSP_EXPORT struct CDoubleArray inputParameters = {params, PARAMS_LENGTH};
const char* params_names[] = {"freq", "q", "vol", "gain"};
DSP_EXPORT struct CStringArray inputParametersNames = {params_names, PARAMS_LENGTH};
//array<string> inputParametersUnits = {};
//double paramMin[] = {0};
//array<double> inputParametersMin(paramMin);
//double paramMax[] = {1};
//array<double> inputParametersMax(paramMax);
//array<double> inputParametersDefault = {};
//array<string> inputParametersFormats = {};
//array<string> inputParametersEnums = {}};
//int steps[] = {-1};
//DSP_EXPORT array<int> inputParametersSteps(steps);

double freq = 0.;
double q = 1.;
double volume = 1.;
double gain = 1.;


DSP_EXPORT bool initialize()
{
    return true;
}

DSP_EXPORT void processBlock(struct BlockData* data)
{
    for(uint c = 0; c < audioInputsCount; ++c)
        for(uint s = 0; s < data->samplesToProcess; ++s)
        {
            double x = data->samples[c][s];
            double y = x;
            data->samples[c][s] = volume*y;
        }
}

DSP_EXPORT void updateInputParameters()
{
    freq = params[0]*params[0];
    q = 4.*params[1] + .707;
    volume = params[2];
    gain = 20.*params[3];
}

/* void updateInputParametersForBlock()
{

} */

DSP_EXPORT void reset()
{
}

/*DSP_EXPORT int getTailSize()
{
	return 0;
} */

/* bool wantSilence()
{
  return false;
} */
