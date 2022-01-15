#include "tfhe/rlwe.h"

void torusPolynomialCopy(TorusPolynomial *result, TorusPolynomial *sample)
{
    const int32_t N = result->N;
    const Torus32 *s = sample->coefsT;
    Torus32 *r = result->coefsT;

    for (int32_t i = 0; i < N; ++i)
        r[i] = s[i];
}

void nonFFTmul(Torus32 *ans, IntPolynomial *S, TorusPolynomial *A, int mod)
{
    for (int i = 0; i < S->N; i++){
        for (int j = 0; j < A->N; j++){
            int64_t temp = ((int64_t)S->coefs[i] * (int64_t)A->coefsT[j]);

            if (i + j < mod)
                ans[i + j] += temp;
            else
                ans[(i + j) % mod] -= temp; 
        }
    }

}

void torusPolynomialSubMulR(TorusPolynomial *result, IntPolynomial *S, TorusPolynomial *A)
{
    int n = result->N;
    for (int i = 0; i < n; i++){
        result->coefsT[i] = -(result->coefsT[i]);
    }

    nonFFTmul(result->coefsT, S, A, n);

    for (int i = 0; i < n; i++){
        result->coefsT[i] = -(result->coefsT[i]);
    }
}

void torusPolynomialAddMulR(TorusPolynomial *result, IntPolynomial *S, TorusPolynomial *A)
{
    int n = result->N;
    
    nonFFTmul(result->coefsT, S, A, n);
}

void tLwePhase(TorusPolynomial *phase, TLweSample *sample, TLweKey *key)
{
    const int32_t k = key->params->k;

    torusPolynomialCopy(phase, sample->b); // phi = b

    for (int32_t i = 0; i < k; ++i)
        torusPolynomialSubMulR(phase, &key->key[i], &sample->a[i]);
}

void tLweApproxPhase(TorusPolynomial *message, TorusPolynomial *phase, int32_t Msize, int32_t N)
{
    for (int32_t i = 0; i < N; ++i)
        message->coefsT[i] = approxPhase(phase->coefsT[i], Msize);
}

void tLweSymDecrypt(TorusPolynomial *result, TLweSample *sample, TLweKey *key, int32_t Msize)
{
    tLwePhase(result, sample, key);
    tLweApproxPhase(result, result, Msize, key->params->N);
}

void torusPolynomialSubTo(TorusPolynomial *result, TorusPolynomial *poly)
{
    int n = result->N;
    for (int i = 0; i < n; i++){
        result->coefsT[i] -= poly->coefsT[i];
    }
}

void torusPolynomialAddTo(TorusPolynomial *result, TorusPolynomial *poly)
{
    int n = result->N;
    for (int i = 0; i < n; i++){
        result->coefsT[i] += poly->coefsT[i];
    }
}