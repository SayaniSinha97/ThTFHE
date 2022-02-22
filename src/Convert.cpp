#include <iostream>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <tfhe/lwe-functions.h>
#include <tfhe/numeric_functions.h>
#include <tfhe/tlwe_functions.h>
#include <random>
#include <time.h>
#include "threshold_decryption_functions.hpp"
#define MSIZE 2

void TLweFromLwe(TLweSample *ring_cipher, LweSample *cipher, TLweParams *tlwe_params){
	int N = tlwe_params->N;
	ring_cipher->a[0].coefsT[0] = cipher->a[0];
	ring_cipher->b->coefsT[0] = cipher->b;
	for(int i = 1; i < N; i++){
		ring_cipher->a[0].coefsT[i] = -cipher->a[N-i];
	}
}

void TLweKeyFromLweKey(const LweKey *lwe_key, TLweKey *tlwe_key){
	int N = tlwe_key->params->N;
	tlwe_key->key[0].coefs[0] = lwe_key->key[0];
	for(int i = 0; i < N; i++){
		tlwe_key->key[0].coefs[i] = lwe_key->key[i];
	}
}

void Evaluate(LweSample *result, LweSample *cipher1, LweSample *cipher2, TFheGateBootstrappingCloudKeySet *cloud_key){
	for(int i = 0; i < 32; i++){
		bootsAND(&result[i], &cipher1[i], &cipher2[i], cloud_key);
	}
}

void BitwiseEncrypt(LweSample *result_array, int msg, TFheGateBootstrappingSecretKeySet *key){
	for (int i = 0; i < 32; i++){
        bootsSymEncrypt(&result_array[i], (msg >> i) & 1, key);
    }
}

int directDecrypt(LweSample *cipher, TFheGateBootstrappingSecretKeySet *key){
	int res = 0;
	for(int i = 0; i < 32; i++){
		res += bootsSymDecrypt(&cipher[i], key) << i;
	}
	return res;
}

int main(){
	FILE *skFile = fopen("test/secret.key", "rb");
    auto key = new_tfheGateBootstrappingSecretKeySet_fromFile(skFile);
    fclose(skFile);

    FILE *plaintext1 = fopen("test/plain22.txt", "r");
    int32_t msg1;
    fscanf(plaintext1, "%d", &msg1);
    fclose(plaintext1);

    FILE *plaintext2 = fopen("test/plain23.txt", "r");
    int32_t msg2;
    fscanf(plaintext2, "%d", &msg2);
    fclose(plaintext2);
    
    FILE *paramFile = fopen("test/secret.params", "rb");
    auto params = new_tfheGateBootstrappingParameterSet_fromFile(paramFile);
    fclose(paramFile);

    FILE *cloudKeyFile = fopen("test/cloud.key", "rb");
    auto cloud_key = new_tfheGateBootstrappingCloudKeySet_fromFile(cloudKeyFile);
    fclose(cloudKeyFile); 

    LweSample *ciphertext1 = new_gate_bootstrapping_ciphertext_array(32, params);
    BitwiseEncrypt(ciphertext1, msg1, key);

    LweSample *ciphertext2 = new_gate_bootstrapping_ciphertext_array(32, params);
    BitwiseEncrypt(ciphertext2, msg2, key);
    LweSample *resultOfEval = new_gate_bootstrapping_ciphertext_array(32, params);
    Evaluate(resultOfEval, ciphertext1, ciphertext2, cloud_key);

    int dmsg = directDecrypt(resultOfEval, key); //Direct decryption of Torus LWE ciphertext

    int32_t resOR = msg1 & msg2;
    std::cout << "\nExpected output: " << resOR << "\n";
    std::cout << "Decrypted value: " << dmsg << "\n";

    TLweParams *tlwe_params = new_TLweParams(1024, 1, 3e-8, 0.2);
    TLweKey *tlwe_key = new_TLweKey(tlwe_params);
	tLweKeyGen(tlwe_key);

	TLweSample *resultOfEvalT;
	TorusPolynomial *result_plaintext;
  TorusPolynomial *direct_result_plaintext;
	TLweKeyFromLweKey(key->lwe_key, tlwe_key);
    shareSecret(3, 5, tlwe_key, tlwe_params);
    double bound = 0.0125;
    int result_msg, direct_result_msg, dbit, tbit;
    std::vector<int> subset{1,2,4};
    while(bound > 1e-5){
	    result_msg = 0;
      direct_result_msg = 0;
	    //Assuming k = 1, n = N = 1024, and converting lwe ciphertext of each of the result bit into corresponding ring-lwe ciphertext one by one and decrypting
		for (int i = 0; i < 32; i++){
			resultOfEvalT = new_TLweSample(tlwe_params);
			TLweFromLwe(resultOfEvalT, &resultOfEval[i], tlwe_params);
			result_plaintext = new_TorusPolynomial(tlwe_params->N);    //result message polynomial after threshold decryption
      direct_result_plaintext = new_TorusPolynomial(tlwe_params->N);    //result message polynomial after direct decryption(for comparison)
      /* Direct decryption of Torus ring-LWE ciphertext*/
      tLwePhase(direct_result_plaintext, resultOfEvalT, tlwe_key);
      dbit = (direct_result_plaintext->coefsT[0] > 0) ? 1 : 0;
      direct_result_msg += dbit << i;
      /* Threshold decryption of Torus ring-LWE ciphertext*/
			thresholdDecrypt(result_plaintext, resultOfEvalT, tlwe_params, subset, 3, 5, bound);
      tbit = result_plaintext->coefsT[0] > 0 ? 1 : 0;
			result_msg += tbit << i;
		}
		std::cout << "\nstandard deviation of Gaussian smudging noise: " << bound << "\n";
    std::cout <<"result after threshold decryption: " << result_msg << "\n";
   std::cout << "result after direct decryption: " << direct_result_msg << "\n";
		bound /= 2;
	}
}