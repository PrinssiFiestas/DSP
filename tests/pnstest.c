// MIT License
// Copyright (c) 2024 Lauri Lorenzo Fiestas
// https://github.com/PrinssiFiestas/DSP/blob/main/LICENCE

#include <dsp/iir.h>

#include <dspapi.h>
#include <chelpers.h>

#include <stdio.h>
#include <stdarg.h>
#include <stdbool.h>

DSP_EXPORT void*          host=null;
DSP_EXPORT HostPrintFunc* hostPrint=null;
DSP_EXPORT double         sampleRate = 0;
DSP_EXPORT uint           audioInputsCount = 0;

void dsp_printf(const char* fmt, ...)
{
    char buf[1024] = "";
    va_list args;
    va_start(args, fmt);
    vsprintf(buf, fmt, args);
    print(buf);
    va_end(args);
}

//-------------------------------------------------------

DSP_EXPORT string const name  = "DSP";
DSP_EXPORT string description = "DSP testing";

enum {
    PARAM_FREQ,
    PARAM_Q,
    PARAM_RIPPLE,
    PARAM_VOL,
    PARAM_GAIN,
    PARAMS_LENGTH
};

double params[PARAMS_LENGTH];
DSP_EXPORT struct CDoubleArray inputParameters = { params, PARAMS_LENGTH };
const char* params_names[PARAMS_LENGTH] = {
    [PARAM_FREQ]   = "freq",
    [PARAM_Q]      = "q",
    [PARAM_RIPPLE] = "ripple",
    [PARAM_VOL]    = "volume",
    [PARAM_GAIN]   = "gain" };
DSP_EXPORT struct CStringArray inputParametersNames = { params_names, PARAMS_LENGTH };
// Uncomment these as needed
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

double freq    = 0.;
double q       = 1.;
double damping = 1.;
double ripple  = 0.;
double volume  = 1.;
double gain    = 1.;

DSP_EXPORT bool initialize(void)
{
    iir_set_sample_rate(sampleRate);
    return true;
}

DSP_EXPORT void processBlock(struct BlockData* data)
{
    for (uint c = 0; c < audioInputsCount; ++c)
        for (uint s = 0; s < data->samplesToProcess; ++s)
        {
            double x = data->samples[c][s];
            double y = x;
            //for (uint i = 0; i < 128; ++i, x = y) // uncomment for performance test
            {
                // TODO TEST EVERYTHING
                static IIRFilter filter[2][IIR_POLES(16)];
                y = iir_fast_high_pass12(filter[c], x, freq, damping);
            }
            data->samples[c][s] = volume*y;
        }
}

DSP_EXPORT void updateInputParameters(void)
{
    freq    = params[PARAM_FREQ]*params[PARAM_FREQ];
    q       = 4.*params[PARAM_Q] + .707;
    ripple  = params[PARAM_RIPPLE];
    damping = 1./q;
    volume  = params[PARAM_VOL];
    gain    = 20.*params[PARAM_GAIN];
}

DSP_EXPORT void reset(void)
{ // TODO reset filter state
}

/*DSP_EXPORT int getTailSize()
{
	return 0;
} */

/* bool wantSilence()
{
  return false;
} */
