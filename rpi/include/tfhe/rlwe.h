#pragma once

#include <stdint.h>
#include "tfhe/common.h"
#include "tfhe/lwe.h"


typedef struct {
   int32_t N;
   Torus32* coefsT;
} TorusPolynomial;

typedef struct {
   TorusPolynomial *a; ///< array of length k+1: mask + right term
   TorusPolynomial *b; ///< alias of a[k] to get the right term
   double current_variance; ///< avg variance of the sample
   int32_t k;
} TLweSample;

typedef struct {
   int32_t N; ///< a power of 2: degree of the polynomials
   int32_t k; ///< number of polynomials in the mask
   double alpha_min; ///< minimal noise s.t. the sample is secure
   double alpha_max; ///< maximal noise s.t. we can decrypt
   LweParams extracted_lweparams; ///< lwe params if one extracts
} TLweParams;

typedef struct {
   int32_t N;
   int32_t* coefs;
} IntPolynomial;

typedef struct {
   TLweParams *params; ///< the parameters of the key
   IntPolynomial *key; ///< the key (i.e k binary polynomials)
} TLweKey;


void tLweSymDecrypt(TorusPolynomial *result, TLweSample *sample, TLweKey *key, int32_t Msize);
void tLweApproxPhase(TorusPolynomial *message, TorusPolynomial *phase, int32_t Msize, int32_t N);
void tLwePhase(TorusPolynomial *phase, TLweSample *sample, TLweKey *key);
void torusPolynomialCopy(TorusPolynomial *result, TorusPolynomial *sample);
void torusPolynomialSubMulR(TorusPolynomial* result, IntPolynomial* poly1, TorusPolynomial* poly2);
void torusPolynomialAddMulR(TorusPolynomial *result, IntPolynomial *poly1, TorusPolynomial *poly2);
void torusPolynomialAddTo(TorusPolynomial *result, TorusPolynomial *poly);
void torusPolynomialSubTo(TorusPolynomial *result, TorusPolynomial *poly);
