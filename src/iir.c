// MIT License
// Copyright (c) 2024 Lauri Lorenzo Fiestas
// https://github.com/PrinssiFiestas/DSP/blob/main/LICENCE

#include <dsp/iir.h>
#include <stdbool.h>
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
    double freq  = normalized_freq * .65; // fix approximation instability
    double alpha = damping*freq*(1. - freq);
    double beta  = (.5 - alpha) / (1. + 2*alpha);
    double cos01 = 2*(freq - 3.)*freq*freq + 1.;
    double gamma = (.5 + beta)*cos01;
    filter[1].a  = .5 + beta - gamma;
    filter[0].a  = filter[1].a/2;
    filter[2].a  = filter[0].a;
    filter[1].b  = 2*gamma;
    filter[2].b  = -2*beta;
}

double fast_low_pass12(Filter filter[3], double input, double normalized_freq, double damping)
{
        coeffs_fast_low_pass12(filter, normalized_freq, damping);
	return apply_filter(filter, 2, input, NULL);
}





typedef struct iir_cheby_coeffs
{
    double a0;
    double a1;
    double a2;
    double b1;
    double b2;
} IIRChebyCoeffs;

static void iir_cheby_coeffs(
    IIRChebyCoeffs* coeff,
    size_t          poles,
    size_t          p,
    double          freq,
    double          ripple,
    bool            is_high_pass)
{
    // Pole location in unit circle
    double rp = -cos(M_PI/(2*poles) + p*M_PI/poles);
    double ip =  sin(M_PI/(2*poles) + p*M_PI/poles);

    // Warp from a circle to an ellipse
    if (ripple != 0.) {
        double r  = 1./(1. - ripple);
        double es = sqrt(r*r - 1.);
        double vx = 1./poles * log(1./es + sqrt(1./(es*es) + 1.));
        double kx = 1./poles * log(1./es + sqrt(1./(es*es) - 1.));
               kx = (exp(kx) + exp(-kx)) / 2;
               rp = rp * ((exp(vx) - exp(-vx)) / 2) / kx;
               ip = ip * ((exp(vx) + exp(-vx)) / 2) / kx;
    }
    // s-domain to z-domain conversion
    double t  = 2*tan(.5);
    double w  = 2*M_PI*freq;
    double m  = rp*rp + ip*ip;
    double d  = 4. - 4*rp*t + m*t*t;
    double x0 = t*t/d;
    double x1 = 2*t*t/d;
    double x2 = t*t/d;
    double y1 = (8. - 2*m*t*t)/d;
    double y2 = (-4. - 4*rp*t - m*t*t)/d;

    // LP to LP, or LP to HP transform
    double k = is_high_pass ?
        -cos(w/2 + .5) / cos(w/2 - .5)
      :  sin(.5 - w/2) / sin(.5 + w/2);
    d = 1. + y1*k - y2*k*k;
    coeff->a0 = (x0 - x1*k + x2*k*k)/d;
    coeff->a1 = (-2*x0*k + x1 + x1*k*k - 2*x2*k)/d;
    coeff->a2 = (x0*k*k - x1*k + x2)/d;
    coeff->b1 = (2*k + y1 + y1*k*k - 2*y2*k)/d;
    coeff->b2 = (-k*k - y1*k + y2)/d;
    if (is_high_pass) {
        coeff->a1 *= -1;
        coeff->b1 *= -1;
    }
}

static void coeffs_chebyshev(
    Filter filter[],
    size_t poles,
    double normalized_freq,
    double normalized_ripple,
    bool   is_high_pass)
{
    double  a[32] = {0};
    double  b[32] = {0};
    double ta[32] = {0};
    double tb[32] = {0};

    double freq   = .5*normalized_freq;
    double ripple = .2925*normalized_ripple;
    a[2] = b[2] = 1.;

    // Loop for each pole pair
    for (size_t p = 0; p < poles/2; ++p) {
        IIRChebyCoeffs coeff = {0};
        iir_cheby_coeffs(&coeff, poles, p, freq, ripple, is_high_pass);

        // Add coefficients to the cascade
        for (size_t i = 0; i <= poles; ++i) {
            ta[i] = a[i];
            tb[i] = b[i];
        }
        for (size_t i = 2; i <= poles + 2; ++i) {
            a[i] = coeff.a0*ta[i] + coeff.a1*ta[i - 1] + coeff.a2*ta[i - 2];
            b[i] =          tb[i] - coeff.b1*tb[i - 1] - coeff.b2*tb[i - 2];
        }
    }
    // Finish combining coefficients
    b[2] = 0.;
    for (size_t i = 0; i <= poles; ++i) {
        a[i] =  a[i + 2];
        b[i] = -b[i + 2];
    }
    // Normalize the gain
    double sa = 0.;
    double sb = 0.;
    for (size_t i = 0; i <= poles; ++i) {
        if ( ! is_high_pass || (i & 1) == 0) {
            sa += a[i];
            sb += b[i];
        } else {
            sa -= a[i];
            sb -= b[i];
        }
    }
    double gain = sa/(1. - sb);
    for (size_t i = 0; i <= poles; ++i) {
        filter[i].a = a[i] / gain;
        filter[i].b = b[i];
    }
}

void coeffs_chebyshev_low_pass(Filter filter[], size_t poles, double freq, double ripple)
{
    // TODO generic freq
    coeffs_chebyshev(filter, poles, freq, ripple, false);
}

void coeffs_chebyshev_high_pass(Filter filter[], size_t poles, double freq, double ripple)
{
    // TODO generic freq
    coeffs_chebyshev(filter, poles, freq, ripple, true);
}

double chebyshev_low_pass(Filter filter[], size_t poles, double input, double freq, double ripple)
{
    coeffs_chebyshev_low_pass(filter, poles, freq, ripple);
    return apply_filter(filter, poles, input, NULL);
}

double chebyshev_high_pass(Filter filter[], size_t poles, double input, double freq, double ripple)
{
    coeffs_chebyshev_high_pass(filter, poles, freq, ripple);
    return apply_filter(filter, poles, input, NULL);
}












static void iir_coeffs_chebyshev(
    IIRCoeffs* filter,
    size_t poles,
    double normalized_freq,
    double normalized_ripple,
    bool   is_high_pass)
{
    double  a[32] = {0};
    double  b[32] = {0};
    double ta[32] = {0};
    double tb[32] = {0};

    double freq   = .5*normalized_freq;
    double ripple = .2925*normalized_ripple;
    a[2] = b[2] = 1.;

    // Loop for each pole pair
    for (size_t p = 0; p < poles/2; ++p) {
        IIRChebyCoeffs coeff = {0};
        iir_cheby_coeffs(&coeff, poles, p, freq, ripple, is_high_pass);

        // Add coefficients to the cascade
        for (size_t i = 0; i <= poles; ++i) {
            ta[i] = a[i];
            tb[i] = b[i];
        }
        for (size_t i = 2; i <= poles + 2; ++i) {
            a[i] = coeff.a0*ta[i] + coeff.a1*ta[i - 1] + coeff.a2*ta[i - 2];
            b[i] =          tb[i] - coeff.b1*tb[i - 1] - coeff.b2*tb[i - 2];
        }
    }
    // Finish combining coefficients
    b[2] = 0.;
    for (size_t i = 0; i <= poles; ++i) {
        a[i] =  a[i + 2];
        b[i] = -b[i + 2];
    }
    // Normalize the gain
    double sa = 0.;
    double sb = 0.;
    for (size_t i = 0; i <= poles; ++i) {
        if ( ! is_high_pass || (i & 1) == 0) {
            sa += a[i];
            sb += b[i];
        } else {
            sa -= a[i];
            sb -= b[i];
        }
    }
    double gain = sa/(1. - sb);
    for (size_t i = 0; i <= poles; ++i) {
        filter->a[i] = a[i] / gain;
        filter->b[i] = b[i];
    }
}

void iir_coeffs_chebyshev_low_pass(IIRCoeffs* filter, size_t poles, double freq, double ripple)
{
    // TODO generic freq
    iir_coeffs_chebyshev(filter, poles, freq, ripple, false);
}

void iir_coeffs_chebyshev_high_pass(IIRCoeffs* filter, size_t poles, double freq, double ripple)
{
    // TODO generic freq
    iir_coeffs_chebyshev(filter, poles, freq, ripple, true);
}

double iir_chebyshev_low_pass(IIRFilter filter[], size_t poles, double input, double freq, double ripple)
{
    iir_coeffs_chebyshev_low_pass(&filter->coeff, poles, freq, ripple);
    return iir_apply_filter(filter, poles, input, NULL);
}

double iir_chebyshev_high_pass(IIRFilter filter[], size_t poles, double input, double freq, double ripple)
{
    iir_coeffs_chebyshev_high_pass(&filter->coeff, poles, freq, ripple);
    return iir_apply_filter(filter, poles, input, NULL);
}


