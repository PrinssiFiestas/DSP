// MIT License
// Copyright (c) 2024 Lauri Lorenzo Fiestas
// https://github.com/PrinssiFiestas/DSP/blob/main/LICENCE

#include <gpc/assert.h>
#include "../src/fft.c"

int main(void)
{
    gp_test("Bit reversal sort"); {
        double res[8]    = { 0, 1, 2, 3, 4, 5, 6, 7 };
        double ims[8]    = { 0, 1, 2, 3, 4, 5, 6, 7 };
        double sorted[8] = { 0, 4, 2, 6, 1, 5, 3, 7 }; // bit patterns reversed
        fft_bit_reversal_sort(res, ims, sizeof res/sizeof res[0]);
        for (size_t i = 0; i < sizeof res/sizeof res[0]; ++i)
            gp_assert(res[i] == sorted[i] && ims[i] == sorted[i], i, res[i], ims[i], sorted[i]);
    }
}
