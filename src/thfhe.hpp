#pragma once

#include <map>
#include <iostream>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <tfhe/lwe-functions.h>
#include <tfhe/numeric_functions.h>
#include <tfhe/tlwe_functions.h>
#include <tfhe/tfhe_garbage_collector.h>
#include <random>
#include <bits/stdc++.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/blas.hpp>
#include <omp.h>
#include <cblas.h>

#define MSIZE 2
#define NSAMPLES 20
namespace ublas = boost::numeric::ublas;
/* Pubkey is just a set of encryptions of 0.
 * b_i = a_i * s + e_i
 * To encrypt:
 * Randomly choose a subset, take sum, add message to sum(b_i).
 */
class ThFHEPubKey {

private:
    std::default_random_engine generator;
    std::uniform_int_distribution<int> *distribution;

public:
    LweSample **samples;
    int n;
    int n_samples;
    double alpha;
    ThFHEPubKey(TFheGateBootstrappingSecretKeySet *sk, int n_samples);

    void Encrypt(LweSample *result, int32_t message);
};

class ThFHEKeyShare {
public:
    std::map<int, TLweKey*> shared_key_repo;	/* Stores <group_id>: <key_share> */
    TLweKey* GetShare(int group_id);
    void PartialDecrypt(TLweSample* ciphertext, TLweParams* params, TorusPolynomial* partial_ciphertext, std::vector<int> parties, int t, int p, double sd);
};

class ThFHE {
private:
    std::map<std::pair<int, int>, int> ncr_cacheT;	            /* Stores <<n, r>: C(n, r)> */
    std::map<std::pair<int, int>, TLweKey*> shared_key_repo;	/* Stores <<party_id, group_id>: key_share> */

    void distributeShares(ublas::matrix<int>& S, int t, int k, int p, TLweParams *params);
    void shareSecret(int t, int p, TLweKey *key, TLweParams *params);


public:
    void KeyGen(int t, int T);
    ThFHEPubKey *pk;
    TFheGateBootstrappingSecretKeySet *sk;
    TFheGateBootstrappingCloudKeySet *bk;    
    void GetShareSet(int party, ThFHEKeyShare *share);    
};

void findParties(std::vector<int>& pt, int gid, int t, int p);
int findGroupId(std::vector<int> parties, int t, int p);
int finalDecrypt(TLweSample* ciphertext, TorusPolynomial** partial_ciphertexts, TLweParams* params, std::vector<int> parties, int t, int p);
void TLweFromLwe(TLweSample *ring_cipher, LweSample *cipher, TLweParams *tlwe_params);
TFheGateBootstrappingParameterSet *initialize_gate_bootstrapping_params();