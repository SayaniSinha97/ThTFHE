#include "tfhe/common.h"
#include "tfhe/lwe.h"

Torus32 lwePhase(const LweSample *sample, const LweKey *key)
{
    const int32_t n = key->params->n;
    Torus32 axs = 0;
    const Torus32 * a = sample->a;
    const int32_t * k = key->key;
    for (int32_t i = 0; i < n; ++i) 
	   axs += a[i]*k[i]; 
    Torus32 ans = sample->b - axs;
    
    return ans;
}

Torus32 approxPhase(Torus32 phase, int32_t Msize){
    uint64_t interv = ((UINT64_C(1)<<63)/Msize)*2;
    uint64_t half_interval = interv/2;
    uint64_t phase64 = (((uint64_t)phase)<<32) + half_interval;
    phase64 -= phase64%interv;
    return (int32_t)(phase64>>32); 
}

int32_t modSwitchFromTorus32(Torus32 phase, int32_t Msize){
    uint64_t interv = ((UINT64_C(1)<<63)/Msize)*2; // width of each intervall
    uint64_t half_interval = interv/2; // begin of the first intervall
    uint64_t phase64 = (((uint64_t)phase)<<32) + half_interval;
    //floor to the nearest multiples of interv
    return phase64/interv;
}


Torus32 lweSymDecrypt(const LweSample *sample, const LweKey *key, const int32_t Msize)
{
    Torus32 phi;
    phi = lwePhase(sample, key);
    return approxPhase(phi, Msize);
}
