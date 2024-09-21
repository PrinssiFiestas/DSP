// MIT License
// Copyright (c) 2024 Lauri Lorenzo Fiestas
// https://github.com/PrinssiFiestas/DSP/blob/main/LICENCE

#include <dsp/iir.h>
#include <math.h>

// Performance notes
//
// Explicit poles parameter instead of encapsulating poles to IIRFilter gives
// the optimizer best chance to unroll loop with different IIR_MAX_POLES values
// in iir_apply_filter().

static double iir_global_sample_rate;
static double iir_global_sample_time;

void iir_set_sample_rate(double sr)
{
    iir_global_sample_rate = sr;
    iir_global_sample_time = 1./sr;
}

double iir_sample_rate(void)
{
    return iir_global_sample_rate;
}

double iir_sample_time(void)
{
    return iir_global_sample_time;
}

double iir_apply_filter(IIRFilter* filter, size_t poles, double input, const IIRNonlinearities* non_lins)
{
    filter->x[0] = input;
    filter->y[0] = filter->coeff.a[0]*input;
    for (size_t i = poles; i > 0; --i) {
        if (non_lins)
            filter->y[0] += filter->coeff.a[i]*non_lins->fx[i](filter->x[i])
                          + filter->coeff.b[i]*non_lins->fy[i](filter->y[i]);
        else
            filter->y[0] += filter->coeff.a[i]*filter->x[i] + filter->coeff.b[i]*filter->y[i];

        filter->x[i] = filter->x[i - 1];
        filter->y[i] = filter->y[i - 1];
    }
    const double anti_denormal = 1e-30;
    filter->y[0] += anti_denormal;
    filter->y[0] -= anti_denormal;

    return filter->y[0];
}

IIRNonlinearities iir_nonlinearities(double(*f)(double))
{
    IIRNonlinearities fs;
    for (size_t i = 0; i < IIR_MAX_POLES; ++i)
        fs.fy[i] = fs.fx[i] = f;
    return fs;
}

//-----------------------------------------------------------------------------
//		FAST LP6

void iir_coeffs_fast_low_pass6(IIRCoeffs* coeff, double normalized_freq)
{
	double a0 = normalized_freq;
	double b1 = 1. - normalized_freq;

	coeff->a[0] = a0;
        coeff->a[1] = 0.;
        coeff->b[0] = 0.;
        coeff->b[1] = b1;
}

double iir_fast_low_pass6(IIRFilter* filter, double input, double normalized_freq)
{
        iir_coeffs_fast_low_pass6(&filter->coeff, normalized_freq);
	return iir_apply_filter(filter, 1, input, NULL);
}

//-----------------------------------------------------------------------------
//		FAST LP12

void iir_coeffs_fast_low_pass12(IIRCoeffs* coeff, double normalized_freq, double damping)
{
    normalized_freq *= .65; // fix approximation instability
    double alpha = damping*normalized_freq*(1. - normalized_freq);
    double beta  = (.5 - alpha) / (1. + 2*alpha);
    double cos01 = 2*(normalized_freq - 3.)*normalized_freq*normalized_freq + 1.;
    double gamma = (.5 + beta)*cos01;
    coeff->a[1] = .5 + beta - gamma;
    coeff->a[0] = coeff->a[1]/2;
    coeff->a[2] = coeff->a[0];
    coeff->b[1] = 2*gamma;
    coeff->b[2] = -2*beta;
}

double iir_fast_low_pass12(IIRFilter* filter, double input, double normalized_freq, double damping)
{
        iir_coeffs_fast_low_pass12(&filter->coeff, normalized_freq, damping);
	return iir_apply_filter(filter, 2, input, NULL);
}










double apply_filter(Filter filter[], size_t poles, double input, const IIRNonlinearities* non_lins)
{
    filter[0].x = input;
    filter[0].y = filter[0].a*input;
    for (size_t i = poles; i > 0; --i) {
        if (non_lins)
            filter[0].y += filter[i].a*non_lins->fx[i](filter[i].x)
                         + filter[i].b*non_lins->fy[i](filter[i].y);
        else
            filter[0].y += filter[i].a*filter[i].x + filter[i].b*filter[i].y;

        filter[i].x = filter[i - 1].x;
        filter[i].y = filter[i - 1].y;
    }
    const double anti_denormal = 1e-30;
    filter[0].y += anti_denormal;
    filter[0].y -= anti_denormal;

    return filter[0].y;
}

//-----------------------------------------------------------------------------
//		LOW PASS 12

void coeffss_low_pass12(Filter filter[3], double freq, double q)
{
	double freq1 = freq <= 1. ? freq : 2.*freq*iir_sample_time();
	double theta = M_PI*freq1;
	double d     = 1./q;
	double dsint = d*sin(theta)/2.;
	double beta  = .5*(1. - dsint)/(1. + dsint);
	double gamma = (.5 + beta)*cos(theta);
	double delta = .5 + beta - gamma;
	filter[0].a = delta/2.;
	filter[1].a = delta;
	filter[2].a = filter[0].a;
	filter[1].b = 2.*gamma;
	filter[2].b = -2.*beta;
}

double low_pass12(Filter filter[3], double input, double freq, double q)
{
	coeffss_low_pass12(filter, freq, q);
	return apply_filter(filter, 2, input, NULL);
}

void coeffs_fast_low_pass6(Filter filter[2], double normalized_freq)
{
	double a0 = normalized_freq;
	double b1 = 1. - normalized_freq;

	filter[0].a = a0;
        filter[0].b = 0.;
        filter[1].a = 0.;
        filter[1].b = b1;
}

double fast_low_pass6(Filter filter[2], double input, double normalized_freq)
{
        coeffs_fast_low_pass6(filter, normalized_freq);
	return apply_filter(filter, 1, input, NULL);
}

//-----------------------------------------------------------------------------
//		FAST LP12

double sin_01(double x)
{
    double x1 = 2*x - 1.;
    return 1. - x1*x1;
}

double cos_01(double x)
{
    double x1 = x - .5;
    return 4*x1*x1*x1 - 3*x1;
}

void coeffs_fast_low_pass12(Filter filter[3], double normalized_freq, double damping)
{
    normalized_freq *= .65; // fix approximation instability
    double alpha = damping*normalized_freq*(1. - normalized_freq);
    double beta  = (.5 - alpha) / (1. + 2*alpha);
    double cos01 = 2*(normalized_freq - 3.)*normalized_freq*normalized_freq + 1.;
    double gamma = (.5 + beta)*cos01;
    filter[1].a = .5 + beta - gamma;
    filter[0].a = filter[1].a/2;
    filter[2].a = filter[0].a;
    filter[1].b = 2*gamma;
    filter[2].b = -2*beta;
}

double fast_low_pass12(Filter filter[3], double input, double normalized_freq, double damping)
{
        coeffs_fast_low_pass12(filter, normalized_freq, damping);
	return apply_filter(filter, 2, input, NULL);
}

