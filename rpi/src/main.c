#include "tfhe/common.h"
#include "tfhe/lwe.h"
#include "tfhe/rlwe.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "utils/perf.h"
#include "utils/random.h"

#define MSIZE 2

#ifndef TP

#define TP 47

#endif

#if TP == 47
#include "tfhe/rlwe_test_data_4_7.h"
int t = 4, p = 7;

#elif TP == 34
#include "tfhe/rlwe_test_data_3_4.h"
int t = 3, p = 4;

#elif TP == 35
#include "tfhe/rlwe_test_data_3_5.h"
int t = 3, p = 5;

#elif TP == 59
#include "tfhe/rlwe_test_data_5_9.h"
int t = 5, p = 9;

#elif TP == 610
#include "tfhe/rlwe_test_data_6_10.h"
int t = 6, p = 10;

#elif TP == 711
#include "tfhe/rlwe_test_data_7_11.h"
int t = 7, p = 11;

#elif TP == 913
#include "tfhe/rlwe_test_data_9_13.h"
int t = 9, p = 13;

#elif TP == 1015
#include "tfhe/rlwe_test_data_10_15.h"
int t = 10, p = 15;

#else

#error "Undefined TP"

#endif



// In this dumbed down version, only the first t parties combine always
// partial_ciphertext must be already initialised to 0s
void partialDecrypt(TLweSample* ciphertext, TLweParams* params, TorusPolynomial* partial_ciphertext, double sd, int party)
{
	if (party > t){
        return;
    }

    int k = params->k;
	int N = params->N;
	int group_id = 0;
	
    IntPolynomial _key = {
        .N = params->N,
        .coefs = __keydata_coef[party][0]
    };

    TLweKey key = {
        .params = params,
        .key = &_key 
    };
    
	for(int j = 0; j < k; j++){
		torusPolynomialAddMulR(partial_ciphertext, key.key, &ciphertext->a[j]);
	}

	for(int j = 0; j < N; j++){
		partial_ciphertext->coefsT[j] += gaussian32(0, sd);
	}
}

int32_t finalDecrypt(TLweSample* ciphertext, TorusPolynomial** partial_ciphertexts, TLweParams* params)
{
	int N = params->N;
	int32_t result_msg = 0;
	TorusPolynomial* result = ciphertext->b;
	for(int i = 0; i < t; i++){
		if(i == 0){
			torusPolynomialSubTo(result, partial_ciphertexts[i]);
		}
		else{
			torusPolynomialAddTo(result, partial_ciphertexts[i]);
		}
	}
	for (int i = 0; i < 32; i++){
		result_msg += (modSwitchFromTorus32(result->coefsT[i], MSIZE) << i);
	}

    return result_msg;
}

int main()
{
    TorusPolynomial _adata = {
        .N = params.N,
        .coefsT = __adata
    };

    Torus32 *__bdatacpy = (Torus32 *)malloc(1024 * sizeof(Torus32));
    for (int i = 0; i < 1024; i++) __bdatacpy[i] = __bdata[i];

    TorusPolynomial _bdata = {
        .N = params.N,
        .coefsT = __bdatacpy
    };

    result.a = &_adata;
    result.b = &_bdata;
    result.k = 1;

    TorusPolynomial plain = {
        .N = 1024,
    };

    plain.coefsT = (Torus32 *)malloc(1024 * sizeof(Torus32));
    for (int i = 0; i < 1024; i++) plain.coefsT[i] = 0;

    printf("Starting decryption for %d out of %d\n", t, p);

    TorusPolynomial **partial_ciphertexts = (TorusPolynomial **)malloc(t * sizeof(TorusPolynomial *));
    for (int i = 0; i < t; i++){
        partial_ciphertexts[i] = (TorusPolynomial *)malloc(sizeof(TorusPolynomial));
        partial_ciphertexts[i]->N = 1024;
        partial_ciphertexts[i]->coefsT = (Torus32 *)malloc(1024 * sizeof(Torus32));
        for (int j = 0; j < 1024; j++) partial_ciphertexts[i]->coefsT[j] = 0;
    }

    double bound = 0.015625;
    for(int i = 0; i < t; i++){
	    clock_t c1 = clock();
	    partialDecrypt(&result, &params, partial_ciphertexts[i], bound, i);
	    clock_t c2 = clock();
	    printf("Partial Decryption Done for party %d Time: %lf ns\n", i, ((double)(c2 - c1)) / CLOCKS_PER_SEC * 1e+9);
    }


    // int32_t c1 = get_cycle_count();
    clock_t c1 = clock();
    int32_t msg = finalDecrypt(&result, partial_ciphertexts, &params);
    clock_t c2 = clock();
    // int32_t c2 = get_cycle_count();
    printf("Message: %d Final Decryption Time: %lf ns\n", msg, ((double)(c2 - c1)) / CLOCKS_PER_SEC * 1e+9);

    return 0;

}
