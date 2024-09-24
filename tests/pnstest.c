// MIT License
// Copyright (c) 2024 Lauri Lorenzo Fiestas
// https://github.com/PrinssiFiestas/DSP/blob/main/LICENCE

#include <dsp/iir.h>

#include <dspapi.h>
#include <chelpers.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>

DSP_EXPORT string const name  = "DSP";
DSP_EXPORT string description = "DSP testing";

DSP_EXPORT void*          host=null;
DSP_EXPORT HostPrintFunc* hostPrint=null;
DSP_EXPORT double         sampleRate=0;
DSP_EXPORT uint           audioInputsCount=0;

// ----------------------------------------------------------------------------
// Debug printing

#if __GNUC__
#define PNS_PRINTF_CHECK_ARGS __attribute__((format(printf,1,2)))
#else
#define PNS_PRINTF_CHECK_ARGS
#endif

void pns_vprintf(const char* fmt, va_list args)
{
    char buf[256]; // Avoid heap usage
    int length = vsnprintf(buf, sizeof buf, fmt, args);
    if (length < 0) {
        print("pns_printf(): output or encoding error.");
        return;
    } else if ((size_t)length > sizeof buf - sizeof"") {
        // This information is not crucial, so it will be only shown if logs are opened.
        print("pns_vprintf(): exceeded max line length.");

        // Give a hint that truncation happened.
        strcpy(buf + sizeof buf - sizeof"...", "...");
    }
    print(buf);
}

PNS_PRINTF_CHECK_ARGS
void pns_printf(const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    pns_vprintf(fmt, args);
    va_end(args);
}

// After calling this with NULL, no printing will be done. This can be used in
// processBlock() to avoid filling logs.
PNS_PRINTF_CHECK_ARGS
void pns_printf_once(const char* fmt, ...)
{
    static bool should_print = true;
    if (fmt == NULL)
        should_print = false;
    if ( ! should_print)
        return;

    va_list args;
    va_start(args, fmt);
    pns_vprintf(fmt, args);
    va_end(args);
}

//-----------------------------------------------------------------------------

typedef enum effect
{
    EFFECT_LOW_PASS6,
    EFFECT_HIGH_PASS6,
    EFFECT_LOW_PASS12,
    EFFECT_HIGH_PASS12,
    EFFECT_BAND_PASS12,
    EFFECT_BAND_STOP12,
    EFFECT_BUTTERWORTH_LOW_PASS12,
    EFFECT_BUTTERWORTH_HIGH_PASS12,
    EFFECT_ALL_PASS6,
    EFFECT_ALL_PASS12,
    EFFECT_LOW_SHELVING,
    EFFECT_HIGH_SHELVING,
    EFFECT_PEAK,
    EFFECT_PEAK_CONST_Q,
    EFFECT_LINKWITZ_RILEY_LOW_PASS12,
    EFFECT_LINKWITZ_RILEY_HIGH_PASS12,
    EFFECT_FAST_LOW_PASS6,
    EFFECT_FAST_HIGH_PASS6,
    EFFECT_FAST_LOW_PASS12,
    EFFECT_FAST_HIGH_PASS12,
    EFFECT_CHEBYSHEV_LOW_PASS,
    EFFECT_CHEBYSHEV_HIGH_PASS,
    EFFECTS_LENGTH,
} Effect;

typedef enum parameter
{
    PARAM_EFFECT,
    PARAM_GAIN,
    PARAM_FREQ,
    PARAM_Q,
    PARAM_POLES,
    PARAM_RIPPLE,
    PARAM_VOL,
    PARAMS_LENGTH
} Parameter;

#define MAX_POLES 16

double params[PARAMS_LENGTH];
const char* params_names[PARAMS_LENGTH] = {
    [PARAM_EFFECT] = "effect",
    [PARAM_FREQ]   = "freq",
    [PARAM_Q]      = "q",
    [PARAM_POLES]  = "poles",
    [PARAM_RIPPLE] = "ripple",
    [PARAM_VOL]    = "volume",
    [PARAM_GAIN]   = "gain",
};
double params_min[PARAMS_LENGTH] = {
    [PARAM_GAIN]   = -18.,
    [PARAM_POLES]  = 2,
};
double params_max[PARAMS_LENGTH] = {
    [PARAM_EFFECT] = EFFECTS_LENGTH - 1 + .01,
    [PARAM_GAIN]   = 18.,
    [PARAM_POLES]  = MAX_POLES,
    [PARAM_RIPPLE] = 30.,
};
const char* params_units[PARAMS_LENGTH] = {
    [PARAM_EFFECT] = "",
    [PARAM_GAIN]   = "dB",
    [PARAM_FREQ]   = "%",
    [PARAM_Q]      = "%",
    [PARAM_POLES]  = "",
    [PARAM_RIPPLE] = "%",
    [PARAM_VOL]    = "%",
};
int params_steps[PARAMS_LENGTH] = {
    [PARAM_EFFECT] = EFFECTS_LENGTH,
    [PARAM_POLES]  = 16 - 2 + 1,
};
const char* params_enums[PARAMS_LENGTH] = { [PARAM_EFFECT] =
    "Low pass 6;"
    "High pass 6;"
    "Low pass 12;"
    "High pass 12;"
    "Band pass 12;"
    "Band stop 12;"
    "Butterworth low pass 12;"
    "Butterworth high pass 12;"
    "All pass 6;"
    "All pass 12;"
    "Low shelving;"
    "High shelving;"
    "Peak;"
    "Peak const q;"
    "Linkwitz Riley low pass 12;"
    "Linkwitz Riley high pass 12;"
    "Fast low pass 6;"
    "Fast high pass 6;"
    "Fast low pass 12;"
    "Fast high pass 12;"
    "Chebyshev low pass;"
    "Chebyshev high pass"
};

DSP_EXPORT struct CDoubleArray inputParameters      = { params,       PARAMS_LENGTH };
DSP_EXPORT struct CStringArray inputParametersNames = { params_names, PARAMS_LENGTH };
DSP_EXPORT struct CDoubleArray inputParametersMin   = { params_min,   PARAMS_LENGTH };
DSP_EXPORT struct CDoubleArray inputParametersMax   = { params_max,   PARAMS_LENGTH };
DSP_EXPORT struct CStringArray inputParametersUnits = { params_units, PARAMS_LENGTH };
DSP_EXPORT struct CIntArray    inputParametersSteps = { params_steps, PARAMS_LENGTH };
DSP_EXPORT struct CStringArray inputParametersEnums = { params_enums, PARAMS_LENGTH };

static Effect effect  = 0;
static double freq    = 0.;
static double q       = 1.;
static double damping = 1.;
static double ripple  = 0.;
static size_t poles   = 2;
static double volume  = 1.;
static double gain    = 1.;
static double gain_db = 0.;

static IIRFilter filter[2/* stereo */][IIR_POLES(MAX_POLES)];

DSP_EXPORT bool initialize(void)
{
    iir_set_sample_rate(sampleRate);
    printf("Hello DSP!");
    return true;
}

static double process_sample(uint ic, double x)
{
    double y = 0.;

    switch (effect) {
    case EFFECT_LOW_PASS6:                  y = iir_low_pass6(filter[ic], x, freq); break;
    case EFFECT_HIGH_PASS6:                 y = iir_high_pass6(filter[ic], x, freq); break;
    case EFFECT_LOW_PASS12:                 y = iir_low_pass12(filter[ic], x, freq, q); break;
    case EFFECT_HIGH_PASS12:                y = iir_high_pass12(filter[ic], x, freq, q); break;
    case EFFECT_BAND_PASS12:                y = iir_band_pass12(filter[ic], x, freq, q); break;
    case EFFECT_BAND_STOP12:                y = iir_band_stop12(filter[ic], x, freq, q); break;
    case EFFECT_BUTTERWORTH_LOW_PASS12:     y = iir_butterworth_low_pass12(filter[ic], x, freq); break;
    case EFFECT_BUTTERWORTH_HIGH_PASS12:    y = iir_butterworth_high_pass12(filter[ic], x, freq); break;
    case EFFECT_ALL_PASS6:                  y = iir_all_pass6(filter[ic], x, freq); break;
    case EFFECT_ALL_PASS12:                 y = iir_all_pass12(filter[ic], x, freq, q); break;
    case EFFECT_LOW_SHELVING:               y = iir_low_shelving(filter[ic], x, freq, gain_db); break;
    case EFFECT_HIGH_SHELVING:              y = iir_high_shelving(filter[ic], x, freq, gain_db); break;
    case EFFECT_PEAK:                       y = iir_peak(filter[ic], x, freq, q, gain_db); break;
    case EFFECT_PEAK_CONST_Q:               y = iir_peak_const_q(filter[ic], x, freq, q, gain_db); break;
    case EFFECT_LINKWITZ_RILEY_LOW_PASS12:  y = iir_linkwitz_riley_low_pass12(filter[ic], x, freq); break;
    case EFFECT_LINKWITZ_RILEY_HIGH_PASS12: y = iir_linkwitz_riley_high_pass12(filter[ic], x, freq); break;
    case EFFECT_FAST_LOW_PASS6:             y = iir_fast_low_pass6(filter[ic], x, freq); break;
    case EFFECT_FAST_HIGH_PASS6:            y = iir_fast_high_pass6(filter[ic], x, freq); break;
    case EFFECT_FAST_LOW_PASS12:            y = iir_fast_low_pass12(filter[ic], x, freq, q); break;
    case EFFECT_FAST_HIGH_PASS12:           y = iir_fast_high_pass12(filter[ic], x, freq, q); break;
    case EFFECT_CHEBYSHEV_LOW_PASS:         y = iir_chebyshev_low_pass(filter[ic], poles, x, freq, ripple); break;
    case EFFECT_CHEBYSHEV_HIGH_PASS:        y = iir_chebyshev_high_pass(filter[ic], poles, x, freq, ripple); break;
    case EFFECTS_LENGTH: break;
    }
    return y;
}

DSP_EXPORT void processBlock(struct BlockData* data)
{
    for (uint ic = 0; ic < audioInputsCount; ++ic) {
        for (uint is = 0; is < data->samplesToProcess; ++is) {
            double y, x = data->samples[ic][is];
            #ifdef PERFORMANCE_TEST
            for (uint i = 0; i < 128; ++i, x = y)
            #endif
                y = process_sample(ic, x);
            data->samples[ic][is] = volume*y;
        }
    }
}

DSP_EXPORT void updateInputParameters(void)
{
    effect  = params[PARAM_EFFECT];
    freq    = params[PARAM_FREQ]*params[PARAM_FREQ];
    q       = 4*params[PARAM_Q] + .707;
    ripple  = params[PARAM_RIPPLE] / 30;
    poles   = params[PARAM_POLES];
    damping = 1./q;
    volume  = params[PARAM_VOL];
    gain_db = params[PARAM_GAIN];
    gain    = 20.*log10(gain_db);
}

DSP_EXPORT void reset(void)
{
    iir_filter_reset(filter[0], MAX_POLES);
    iir_filter_reset(filter[1], MAX_POLES);
}

/*DSP_EXPORT int getTailSize()
{
	return 0;
} */

/* bool wantSilence()
{
  return false;
} */
