// MIT License
// Copyright (c) 2024 Lauri Lorenzo Fiestas
// https://github.com/PrinssiFiestas/DSP/blob/main/LICENCE

#include <dsp/iir.h>
#include <stdbool.h>
#include <math.h>

#define IIR_E           2.7182818284590452354
#define IIR_LOG2E       1.4426950408889634074
#define IIR_LOG10E      0.43429448190325182765
#define IIR_LN2         0.69314718055994530942
#define IIR_LN10        2.30258509299404568402
#define IIR_PI          3.14159265358979323846
#define IIR_PI_2        1.57079632679489661923
#define IIR_PI_4        0.78539816339744830962
#define IIR_1_PI        0.31830988618379067154
#define IIR_2_PI        0.63661977236758134308
#define IIR_2_SQRTPI    1.12837916709551257390
#define IIR_SQRT2       1.41421356237309504880
#define IIR_SQRT1_2     0.70710678118654752440

double iir_sample_rate;

void iir_filter_reset(IIRFilter filter[], size_t poles)
{
    for (size_t i = 0; i <= poles; ++i)
        filter[i].x = filter[i].y = 0.;
}

double iir_apply_filter(IIRFilter filter[], size_t poles, double input, const IIRNonlinearities non_lins[])
{
    filter[0].x = input;
    filter[0].y = filter[0].a*input;
    for (size_t i = poles; i > 0; --i) {
        if (non_lins)
            filter[0].y += filter[i].a*non_lins[i].fx(filter[i].x)
                         + filter[i].b*non_lins[i].fy(filter[i].y);
        else
            filter[0].y += filter[i].a*filter[i].x + filter[i].b*filter[i].y;

        filter[i].x = filter[i - 1].x;
        filter[i].y = filter[i - 1].y;
    }
    const double   anti_denormal = 1e-30;
    filter[0].y += anti_denormal;
    filter[0].y -= anti_denormal;

    return filter[0].y;
}

// ----------------------------------------------------------------------------
// Low Pass 6

void iir_coeffs_low_pass6(IIRFilter filter[IIR_POLES(1)], double freq)
{
    double freq1 = freq <= 1. ? freq : 2*freq/iir_sample_rate;
    double gamma = 2. - cos(IIR_PI*freq1);
    filter[1].b  = gamma - sqrt(gamma*gamma - 1.);
    filter[0].a  = 1. - filter[1].b;
    filter[1].a  = 0.;
}

double iir_low_pass6(IIRFilter filter[IIR_POLES(1)], double input, double freq)
{
    iir_coeffs_low_pass6(filter, freq);
    return iir_apply_filter(filter, 1, input, NULL);
}

// ----------------------------------------------------------------------------
// High Pass 6

void iir_coeffs_high_pass6(IIRFilter filter[IIR_POLES(1)], double freq)
{
    double freq1 = freq <= 1. ? freq : 2*freq/iir_sample_rate;
    double theta = IIR_PI*freq1;
    double gamma = cos(theta)/(1. + sin(theta));
    filter[0].a  = (1. + gamma)/2;
    filter[1].a  = -filter[0].a;
    filter[1].b  = gamma;
}

double iir_high_pass6(IIRFilter filter[IIR_POLES(1)], double input, double freq)
{
    iir_coeffs_high_pass6(filter, freq);
    return iir_apply_filter(filter, 1, input, NULL);
}

//-----------------------------------------------------------------------------
// Low Pass 12

void iir_coeffs_low_pass12(IIRFilter filter[IIR_POLES(2)], double freq, double q)
{
    double freq1 = freq <= 1. ? freq : 2*freq/iir_sample_rate;
    double theta = IIR_PI*freq1;
    double d     = 1./q;
    double dsint = d*sin(theta)/2;
    double beta  = .5*(1. - dsint)/(1. + dsint);
    double gamma = (.5 + beta)*cos(theta);
    double delta = .5 + beta - gamma;
    filter[0].a  = delta/2;
    filter[1].a  = delta;
    filter[2].a  = filter[0].a;
    filter[1].b  = 2*gamma;
    filter[2].b  = -2*beta;
}

double iir_low_pass12(IIRFilter filter[IIR_POLES(2)], double input, double freq, double q)
{
    iir_coeffs_low_pass12(filter, freq, q);
    return iir_apply_filter(filter, 2, input, NULL);
}

// ----------------------------------------------------------------------------
// High Pass 12

void iir_coeffs_high_pass12(IIRFilter filter[IIR_POLES(2)], double freq, double q)
{
    double freq1 = freq <= 1. ? freq : 2*freq/iir_sample_rate;
    double theta = IIR_PI*freq1;
    double d     = 1./q;
    double dsint = d*sin(theta)/2.;
    double beta  = .5*(1. - dsint)/(1. + dsint);
    double gamma = (.5 + beta)*cos(theta);
    double sigma = .5 + beta + gamma;
    filter[0].a  = sigma/2;
    filter[1].a  = -sigma;
    filter[2].a  = filter[0].a;
    filter[1].b  = 2*gamma;
    filter[2].b  = -2*beta;
}

double iir_high_pass12(IIRFilter filter[IIR_POLES(2)], double input, double freq, double q)
{
    iir_coeffs_high_pass12(filter, freq, q);
    return iir_apply_filter(filter, 2, input, NULL);
}

// ----------------------------------------------------------------------------
// Band Pass 12

void iir_coeffs_band_pass12(IIRFilter filter[IIR_POLES(2)], double freq, double q)
{
    double freq1 = freq <= 1. ? freq/2 : freq/iir_sample_rate;
    double k     = tan(IIR_PI*freq1);
    double delta = k*k*q + k + q;
    filter[0].a  = k/delta;
    filter[1].a  = 0.;
    filter[2].a  = -filter[0].a;
    filter[1].b  = -2*q*(k*k - 1.)/delta;
    filter[2].b  = -(k*k*q - k + q)/delta;
}

double iir_band_pass12(IIRFilter filter[IIR_POLES(2)], double input, double freq, double q)
{
    iir_coeffs_band_pass12(filter, freq, q);
    return iir_apply_filter(filter, 2, input, NULL);
}

// ----------------------------------------------------------------------------
// Band Stop 12

void iir_coeffs_band_stop12(IIRFilter filter[IIR_POLES(2)], double freq, double q)
{
    double freq1 = freq <= 1. ? freq/2 : freq/iir_sample_rate;
    double k     = tan(IIR_PI*freq1);
    double delta = k*k*q + k + q;
    filter[0].a  = q*(k*k + 1.)/delta;
    filter[1].a  = 2*q*(k*k - 1.)/delta;
    filter[2].a  = filter[0].a;
    filter[1].b  = -filter[1].a;
    filter[2].b  = -(k*k*q - k + q)/delta;
}

double iir_band_stop12(IIRFilter filter[IIR_POLES(2)], double input, double freq, double q)
{
    iir_coeffs_band_stop12(filter, freq, q);
    return iir_apply_filter(filter, 2, input, NULL);
}

// ----------------------------------------------------------------------------
// Butterworth Low Pass 12

void iir_coeffs_butterworth_low_pass12(IIRFilter filter[IIR_POLES(2)], double freq)
{
    double freq1   = freq <= 1. ? freq/2 : freq/iir_sample_rate;
    double bounded = fmin(fmax(freq1, 1e-5), .5 - 1e-4); // prevent singularities
    double c       = 1./tan(IIR_PI*bounded);
    filter[0].a    = 1./(1. + IIR_SQRT2*c + c*c);
    filter[1].a    = 2*filter[0].a;
    filter[2].a    = filter[0].a;
    filter[1].b    = -2*filter[0].a*(1. - c*c);
    filter[2].b    = -filter[0].a*(1. - IIR_SQRT2*c + c*c);
}

double iir_butterworth_low_pass12(IIRFilter filter[IIR_POLES(2)], double input, double freq)
{
    iir_coeffs_butterworth_low_pass12(filter, freq);
    return iir_apply_filter(filter, 2, input, NULL);
}

// ----------------------------------------------------------------------------
// Butterworth High Pass 12

void iir_coeffs_butterworth_high_pass12(IIRFilter filter[IIR_POLES(2)], double freq)
{
    double freq1   = freq <= 1. ? freq/2 : freq/iir_sample_rate;
    double _freq1  = .5 - freq1;
    double bounded = fmin(fmax(_freq1, 1e-5), .5 - 1e-4); // prevent singularities
    double c       = 1./tan(IIR_PI*bounded);
    filter[0].a    = 1./(1. + IIR_SQRT2*c + c*c);
    filter[1].a    = -2*filter[0].a;
    filter[2].a    = filter[0].a;
    filter[1].b    = -2*filter[0].a*(c*c - 1.);
    filter[2].b    = -filter[0].a*(1. - IIR_SQRT2*c + c*c);
}

double iir_butterworth_high_pass12(IIRFilter filter[IIR_POLES(2)], double input, double freq)
{
    iir_coeffs_butterworth_high_pass12(filter, freq);
    return iir_apply_filter(filter, 2, input, NULL);
}

// ----------------------------------------------------------------------------
// All Pass 6

void iir_coeffs_all_pass6(IIRFilter filter[IIR_POLES(1)], double freq)
{
    double freq1 = freq <= 1. ? freq/2 : freq/iir_sample_rate;
    double theta = IIR_PI*freq1;
    double alpha = (tan(theta) - 1.)/(tan(theta) + 1.);
    filter[0].a  = alpha;
    filter[1].a  = 1.;
    filter[1].b  = -alpha;
}

double iir_all_pass6(IIRFilter filter[IIR_POLES(1)], double input, double freq)
{
    iir_coeffs_all_pass6(filter, freq);
    return iir_apply_filter(filter, 1, input, NULL);
}

// ----------------------------------------------------------------------------
// All Pass 12

void iir_coeffs_all_pass12(IIRFilter filter[IIR_POLES(2)], double freq, double q)
{
    double freq1 = freq <= 1. ? freq/2 : freq/iir_sample_rate;
    double theta = IIR_PI*freq1;
    double alpha = (tan(theta/q) - 1.)/(tan(theta/q) + 1.);
    double beta  = -cos(2*theta);
    filter[0].a  = -alpha;
    filter[1].a  = beta*(1. - alpha);
    filter[2].a  = 1.;
    filter[1].b  = -filter[1].a;
    filter[2].b  = alpha;
}

double iir_all_pass12(IIRFilter filter[IIR_POLES(2)], double input, double freq, double q)
{
    iir_coeffs_all_pass12(filter, freq, q);
    return iir_apply_filter(filter, 2, input, NULL);
}

// ----------------------------------------------------------------------------
// Low Shelving

void iir_coeffs_low_shelving(IIRFilter filter[IIR_POLES(1)], double freq, double gain)
{
    double freq1 = freq <= 1. ? freq/2 : freq/iir_sample_rate;
    double theta = IIR_PI*freq1;
    double beta  = 4./(1. + gain);
    double delta = beta*tan(theta);
    double gamma = (1. - delta)/(1. + delta);
    filter[0].a  = (1. - gamma)/2;
    filter[1].a  = filter[0].a;
    filter[1].b  = gamma;
}

double iir_low_shelving(IIRFilter filter[IIR_POLES(1)], double input, double freq, double gain_db)
{
    double gain = pow(10., gain_db/20);
    iir_coeffs_low_shelving(filter, freq, gain);
    return input + (gain - 1.)*iir_apply_filter(filter, 1, input, NULL);
}

// ----------------------------------------------------------------------------
// High Shelving

void iir_coeffs_high_shelving(IIRFilter filter[IIR_POLES(1)], double freq, double gain)
{
    double freq1 = freq <= 1. ? freq/2 : freq/iir_sample_rate;
    double theta = IIR_PI*freq1;
    double beta  = (1. + gain)/4;
    double delta = beta*tan(theta);
    double gamma = (1. - delta)/(1. + delta);
    filter[0].a  = (1. + gamma)/2;
    filter[1].a  = -filter[0].a;
    filter[1].b  = gamma;
}

double iir_high_shelving(IIRFilter filter[IIR_POLES(1)], double input, double freq, double gain_db)
{
    double gain = pow(10., gain_db/20);
    iir_coeffs_high_shelving(filter, freq, gain);
    return input + (gain - 1.)*iir_apply_filter(filter, 1, input, NULL);
}

// ----------------------------------------------------------------------------
// Peak

void iir_coeffs_peak(IIRFilter filter[IIR_POLES(2)], double freq, double q, double gain)
{
    double freq1 = freq <= 1. ? freq/2 : freq/iir_sample_rate;
    double theta = 2*IIR_PI*freq1;
    double sigma = 4./(1. + gain);
    double beta  = .5*(1. - sigma*tan(.5*theta/q))/(1. + sigma*tan(.5*theta/q));
    double gamma = (.5 + beta)*cos(theta);
    filter[0].a  = .5 - beta;
    filter[1].a  = 0.;
    filter[2].a  = -filter[0].a;
    filter[1].b  = 2*gamma;
    filter[2].b  = -2*beta;
}

double iir_peak(IIRFilter filter[IIR_POLES(2)], double input, double freq, double q, double gain_db)
{
    double gain = pow(10., gain_db/20);
    iir_coeffs_peak(filter, freq, gain, q);
    return input + (gain - 1.)*iir_apply_filter(filter, 2, input, NULL);
}

// ----------------------------------------------------------------------------
// Peak Const Q

void iir_coeffs_peak_const_q(IIRFilter filter[IIR_POLES(2)], double freq, double q, double gain)
{
    double freq1 = freq <= 1. ? freq/2 : freq/iir_sample_rate;
    double theta = IIR_PI*freq1;
    double k     = tan(theta);
    double d0    = 1. + k/q + k*k;
    double e0    = 1. + k/(gain*q) + k*k;
    double alpha = 1. + gain*k/q + k*k;
    double beta  = 2*(k*k - 1.);
    double gamma = 1. - gain*k/q + k*k;
    double delta = 1. - k/q + k*k;
    double eta   = 1. - k/(gain*q) + k*k;

    int boost   = gain >= 1.;
    filter[0].a = boost ? alpha/d0  : d0/e0;
    filter[1].a = boost ? beta/d0   : beta/e0;
    filter[2].a = boost ? gamma/d0  : delta/e0;
    filter[1].b = boost ? -beta/d0  : -beta/e0;
    filter[2].b = boost ? -delta/d0 : -eta/e0;
}

double iir_peak_const_q(IIRFilter filter[IIR_POLES(2)], double input, double freq, double q, double gain_db)
{
    double gain = pow(10., gain_db/20);
    iir_coeffs_peak_const_q(filter, freq, q, gain);
    return iir_apply_filter(filter, 2, input, NULL);
}

// ----------------------------------------------------------------------------
// Linkwitz Riley Low Pass 12

void iir_coeffs_linkwitz_riley_low_pass12(IIRFilter filter[IIR_POLES(2)], double freq)
{
    double sample_rate = freq <= 1. ? 2. : iir_sample_rate;
    double theta = IIR_PI*freq/sample_rate;
    double omega = IIR_PI*freq;
    double kappa = omega/tan(theta);
    double delta = kappa*kappa + omega*omega + 2*kappa*omega;
    filter[0].a  = omega*omega/delta;
    filter[1].a  = 2*filter[0].a;
    filter[2].a  = filter[0].a;
    filter[1].b  = (2*kappa*kappa - 2*omega*omega)/delta;
    filter[2].b  = (2*kappa*omega - kappa*kappa - omega*omega)/delta;
}

double iir_linkwitz_riley_low_pass12(IIRFilter filter[IIR_POLES(2)], double input, double freq)
{
    iir_coeffs_linkwitz_riley_low_pass12(filter, freq);
    return iir_apply_filter(filter, 2, input, NULL);
}

// ----------------------------------------------------------------------------
// Linkwitz Riley High Pass 12

void iir_coeffs_linkwitz_riley_high_pass12(IIRFilter filter[IIR_POLES(2)], double freq)
{
    double sample_rate = freq <= 1. ? 2. : iir_sample_rate;
    double theta = IIR_PI*freq/sample_rate;
    double omega = IIR_PI*freq;
    double kappa = omega/tan(theta);
    double delta = kappa*kappa + omega*omega + 2*kappa*omega;
    filter[0].a  = kappa*kappa/delta;
    filter[1].a  = -2*filter[0].a;
    filter[2].a  = filter[0].a;
    filter[1].b  = (2*kappa*kappa - 2*omega*omega)/delta;
    filter[2].b  = (2*kappa*omega - kappa*kappa - omega*omega)/delta;
}

double iir_linkwitz_riley_high_pass12(IIRFilter filter[IIR_POLES(2)], double input, double freq)
{
    iir_coeffs_linkwitz_riley_high_pass12(filter, freq);
    return iir_apply_filter(filter, 2, input, NULL);
}

// ----------------------------------------------------------------------------
// Fast Low Pass 6

void iir_coeffs_fast_low_pass6(IIRFilter filter[2], double normalized_freq)
{
    double a0 = normalized_freq;
    double b1 = 1. - normalized_freq;

    filter[0].a = a0;
    filter[1].a = 0.;
    filter[1].b = b1;
}

double iir_fast_low_pass6(IIRFilter filter[2], double input, double normalized_freq)
{
    iir_coeffs_fast_low_pass6(filter, normalized_freq);
    return iir_apply_filter(filter, 1, input, NULL);
}

// ----------------------------------------------------------------------------
// Fast High Pass 6

void iir_coeffs_fast_high_pass6(IIRFilter filter[IIR_POLES(1)], double normalized_freq)
{
    double gamma = 1. - 2*normalized_freq;
    filter[0].a  = (1. + gamma)/2;
    filter[1].a  = -filter[0].a;
    filter[1].b  = gamma;
}

double iir_fast_high_pass6(IIRFilter filter[IIR_POLES(1)], double input, double normalized_freq)
{
    iir_coeffs_fast_high_pass6(filter, normalized_freq);
    return iir_apply_filter(filter, 1, input, NULL);
}

//-----------------------------------------------------------------------------
// Fast Low Pass 12

void iir_coeffs_fast_low_pass12(IIRFilter filter[3], double normalized_freq, double damping)
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

double iir_fast_low_pass12(IIRFilter filter[3], double input, double normalized_freq, double damping)
{
    iir_coeffs_fast_low_pass12(filter, normalized_freq, damping);
    return iir_apply_filter(filter, 2, input, NULL);
}

// ----------------------------------------------------------------------------
// Fast High Pass 12

void iir_coeffs_fast_high_pass12(IIRFilter filter[IIR_POLES(2)], double normalized_freq, double damping)
{
    double freq  = normalized_freq * .65; // fix approximation instability
    double alpha = damping*freq*(1. - freq);
    double beta  = (.5 - alpha) / (1. + 2*alpha);
    double cos01 = 2*(freq - 3.)*freq*freq + 1.;
    double gamma = (.5 + beta)*cos01;
    filter[1].a  = .5 + beta + gamma;
    filter[0].a  = filter[1].a/2;
    filter[2].a  = filter[0].a;
    filter[1].a *= -1.;
    filter[1].b  = 2*gamma;
    filter[2].b  = -2*beta;
}

double iir_fast_high_pass12(IIRFilter filter[3], double input, double normalized_freq, double damping)
{
    iir_coeffs_fast_high_pass12(filter, normalized_freq, damping);
    return iir_apply_filter(filter, 2, input, NULL);
}

// -------------------------------------------------------------------.--------
// Chebyshev

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
    double rp = -cos(IIR_PI/(2*poles) + p*IIR_PI/poles);
    double ip =  sin(IIR_PI/(2*poles) + p*IIR_PI/poles);

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
    double w  = 2*IIR_PI*freq;
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

static void iir_coeffs_chebyshev(
    IIRFilter filter[],
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

void iir_coeffs_chebyshev_low_pass(IIRFilter filter[], size_t poles, double freq, double ripple)
{
    double freq1 = freq <= 1. ? freq : 2.*freq/iir_sample_rate;
    iir_coeffs_chebyshev(filter, poles, freq1, ripple, false);
}

void iir_coeffs_chebyshev_high_pass(IIRFilter filter[], size_t poles, double freq, double ripple)
{
    double freq1 = freq <= 1. ? freq : 2.*freq/iir_sample_rate;
    iir_coeffs_chebyshev(filter, poles, freq1, ripple, true);
}

double iir_chebyshev_low_pass(IIRFilter filter[], size_t poles, double input, double freq, double ripple)
{
    iir_coeffs_chebyshev_low_pass(filter, poles, freq, ripple);
    return iir_apply_filter(filter, poles, input, NULL);
}

double iir_chebyshev_high_pass(IIRFilter filter[], size_t poles, double input, double freq, double ripple)
{
    iir_coeffs_chebyshev_high_pass(filter, poles, freq, ripple);
    return iir_apply_filter(filter, poles, input, NULL);
}
