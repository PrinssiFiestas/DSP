// MIT License
// Copyright (c) 2024 Lauri Lorenzo Fiestas
// https://github.com/PrinssiFiestas/DSP/blob/main/LICENCE

#include <dsp/iir.h>
#include <stddef.h>

double iir_apply_filter(IIRFilter* filter, double input, IIRNonlinearities* non_lins)
{
    filter->x[0] = input;
    filter->y[0] = filter->coeff.a[0] * input;
    for (size_t i = IIR_MAX_POLES; i > 0; --i) {
        if (non_lins)
            filter->y[0] += filter->coeff.a[i] * non_lins->fx[i](filter->x[i])
                          + filter->coeff.b[i] * non_lins->fy[i](filter->y[i]);
        else
            filter->y[0] += filter->coeff.a[i] * filter->x[i] + filter->coeff.b[i] * filter->y[i];

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

