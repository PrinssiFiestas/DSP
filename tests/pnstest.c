// MIT License
// Copyright (c) 2024 Lauri Lorenzo Fiestas
// https://github.com/PrinssiFiestas/DSP/blob/main/LICENCE

#include <dsp/iir.h>

#include <dspapi.h>
#include <chelpers.h>

#include <stdbool.h>

DSP_EXPORT void*          host=null;
DSP_EXPORT HostPrintFunc* hostPrint=null;
DSP_EXPORT double         sampleRate = 0;
DSP_EXPORT uint           audioInputsCount = 0;

//-------------------------------------------------------

DSP_EXPORT string const name  = "DSP";
DSP_EXPORT string description = "DSP testing";

enum {
    PARAM_FREQ,
    PARAM_Q,
    PARAM_VOL,
    PARAM_GAIN,
    PARAMS_LENGTH
};

double params[PARAMS_LENGTH];
DSP_EXPORT struct CDoubleArray inputParameters = {params, PARAMS_LENGTH};
const char* params_names[PARAMS_LENGTH] = {"freq", "q", "vol", "gain"};
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
double damping = 1.;
double volume = 1.;
double gain = 1.;

DSP_EXPORT bool initialize(void)
{
    iir_set_sample_rate(sampleRate);
    return true;
}

DSP_EXPORT void processBlock(struct BlockData* data)
{
    static IIRFilter filter[2] = {0};

    for (uint c = 0; c < audioInputsCount; ++c)
        for (uint s = 0; s < data->samplesToProcess; ++s)
        {
            double x = data->samples[c][s];
            double y = x;
            //for (uint i = 0; i < 1024; ++i, x = y) // uncomment for performance test
            {
                #if 1
                //y = iir_fast_low_pass12(&filter[c], x, freq, q);
                y = iir_fast_low_pass6(&filter[c], x, freq);
                #else
                (void)filter;
                static Filter fltr[2][IIR_POLES(2)];
                y = fast_low_pass6(fltr[c], x, freq);
                //y = fast_low_pass12(fltr[c], x, freq, damping);
                //y = low_pass12(fltr[c], x, freq, q);
                #endif
            }
            data->samples[c][s] = volume*y;
        }
}

void coeffs_fast_low_pass12(Filter filter[3], double, double);

DSP_EXPORT void updateInputParameters(void)
{
    freq = params[0]*params[0];
    q = 4.*params[1] + .707;
    damping = 1./q;
    volume = params[2];
    gain = 20.*params[3];
}

/* void updateInputParametersForBlock()
{

} */

DSP_EXPORT void reset(void)
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
