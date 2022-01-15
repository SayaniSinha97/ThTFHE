#pragma once

#include "tfhe/common.h"

typedef struct
{
    const int32_t n;
    const double alpha_min;
    const double alpha_max;
} LweParams;

typedef struct
{
    Torus32 *a;
    Torus32 b;
    double current_variance;
} LweSample;

typedef struct
{
    const LweParams *params;
    int32_t *key;
} LweKey;


// TODO: Port IO functions

Torus32 lweSymDecrypt(const LweSample *sample, const LweKey *key, const int32_t Msize);
Torus32 lwePhase(const LweSample *sample, const LweKey *key);
Torus32 approxPhase(Torus32 phase, int32_t Msize);
int32_t modSwitchFromTorus32(Torus32 phase, int32_t Msize);