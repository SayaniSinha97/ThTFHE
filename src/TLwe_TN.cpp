#include <iostream>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <tfhe/lwe-functions.h>
#include <tfhe/numeric_functions.h>
#include <tfhe/tlwe_functions.h>
#include <random>
#include <time.h>
#include <cstdint>
#include <x86intrin.h>
#include <bits/stdc++.h>
#include "threshold_decryption_functions.hpp"


#define MSIZE 2

int main(int argc, char *argv[]){
	if(argc < 4){
		std::cout << "Please provide values of t, p and party-ids of collaborating t parties as space separated integers in the command line for t-out-of-p threshold decryption.\n";
		return 0;
	}
	int t = atoi(argv[1]);
	int p = atoi(argv[2]);
	int party_id;
	std::vector<int> subset;
	for(int i = 3; i < argc; i++){
		party_id = atoi(argv[i]);
		if(party_id <= p)		/* Otherwise, the party id is invalid */
			subset.push_back(atoi(argv[i]));
	}

	/* Check uniqueness and correctness of the provided party-ids */
	std::sort(subset.begin(), subset.end());
	std::vector<int>::iterator it;
	it = std::unique(subset.begin(), subset.end());
	subset.resize(std::distance(subset.begin(), it));
	if(subset.size() < t){
		std::cout << "Please provide at least " << t << " correct and unique party-ids to get result of " << t << "-out-of-" << p << " threshold decrypton.\n";
		return 0;
	}

	/* Read from plaintext file */
	FILE *plaintext = fopen("test/plain22.txt", "r");
    int32_t msg;
    fscanf(plaintext, "%d", &msg);
    fclose(plaintext);
    std::cout << "Plaintext: " << msg <<"\n";

    /* Set Up */
	TLweParams *params = new_TLweParams(1024, 2, 0.01, 0.2);
	TLweKey *key = new_TLweKey(params);
	tLweKeyGen(key);

	/* Encryption */
	TLweSample *ciphertext = new_TLweSample(params);
	TorusPolynomial* mu = new_TorusPolynomial(params->N);
	for(int i = 0; i < params->N; i++){
		mu->coefsT[i] = 0;
	}
	for(int i = 0; i < 32; i++){
		mu->coefsT[i] += modSwitchToTorus32((msg >> i) & 1, MSIZE);
	}
	tLweSymEncrypt(ciphertext, mu, 3e-8, key);

    /* Direct Decryption */
    int dmsg = 0;
    TorusPolynomial* res = new_TorusPolynomial(params->N);
    tLweSymDecrypt(res, ciphertext, key, MSIZE);
    for(int i = 0; i < 32; i++){
		dmsg += (modSwitchFromTorus32(res->coefsT[i], MSIZE) << i);
    }
    std::cout << "Direct Decryption result: " << dmsg << std::endl;

    /* Threshold Decryption */

	// struct timespec share_start = {0, 0};
	// struct timespec share_end = {0, 0};
	// clock_gettime(CLOCK_MONOTONIC, &share_start);	/*measure time in seconds*/
	unsigned int high, low;
	__asm__ __volatile__("xorl %%eax,%%eax\n cpuid \n" ::: "%eax", "%ebx", "%ecx", "%edx");
    __asm__ __volatile__ ("rdtsc" : "=a" (low), "=d" (high));
	auto clock_start_sharing = (static_cast<uint64_t>(high) << 32) | low;	/*measure time in clock cycles*/

    shareSecret2(t, p, key, params);

    __asm__ __volatile__ ("rdtsc" : "=a" (low), "=d" (high));
    auto clock_stop_sharing = (static_cast<uint64_t>(high) << 32) | low;
	// clock_gettime(CLOCK_MONOTONIC, &share_end);

	// std::cout << "Secret sharing time:" << (((double)share_end.tv_nsec + 1.0e+9 * share_end.tv_sec) - ((double)share_start.tv_nsec + 1.0e+9 * share_start.tv_sec)) * 1.0e-9 << "sec" << std::endl;
    auto cycle_count_sharing = clock_stop_sharing - clock_start_sharing;
    std::cout << "Secret Sharing Time: " << cycle_count_sharing << " cycles" << std::endl;
    // struct timespec start_time = {0, 0};
    // struct timespec end_time = {0, 0};
    
    int rbit;
    int result_msg;
    TorusPolynomial* result_plaintext = new_TorusPolynomial(params->N);
    double bound = 0.05;
    unsigned int lo,hi;

 //    while(bound > 1e-3){
	//     result_msg = 0;
	//     // clock_gettime(CLOCK_MONOTONIC, &start_time);
	//     __asm__ __volatile__("xorl %%eax,%%eax\n cpuid \n" ::: "%eax", "%ebx", "%ecx", "%edx");
	//     __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	// 	auto clock_start_decryption = (static_cast<uint64_t>(hi) << 32) | lo;
	// 	thresholdDecrypt(result_plaintext, ciphertext, params, subset, t, p, bound);
	// 	for (int i = 0; i < 32; i++){
	// 		result_msg += (modSwitchFromTorus32(result_plaintext->coefsT[i], MSIZE) << i);
	// 	}
		
	// 	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	// 	auto clock_stop_decryption = (static_cast<uint64_t>(hi) << 32) | lo;
	//     // clock_gettime(CLOCK_MONOTONIC, &end_time);
	//     auto cycle_count_decryption = clock_stop_decryption - clock_start_decryption;
	//     // std::cout << t << "-out-of-" << p << " Threshold Decryption result(with bound " << bound << "): "<< result_msg << ". Decryption Time: " << (((double)end_time.tv_nsec + 1.0e+9 * end_time.tv_sec) - ((double)start_time.tv_nsec + 1.0e+9 * start_time.tv_sec)) * 1.0e-9 << " sec    " << std::endl;
	//     std::cout << t << "-out-of-" << p << " Threshold Decryption result(with bound " << bound << "): "<< result_msg << ". Decryption Time: " << cycle_count_decryption << " cycles" << std::endl;
	//     bound /= 2;
	// }
	uint64_t* cycle_counts_partial = new uint64_t[t];
	uint64_t* cycle_counts_final = new uint64_t[t];
	TorusPolynomial** partial_ciphertexts = new TorusPolynomial*[t];
	for(int i = 0; i < t; i++){
		partial_ciphertexts[i] = new_TorusPolynomial(params->N);
	}
	// TorusPolynomial** result_plaintexts = new TorusPolynomial*[t];
	// for(int i = 0; i < t; i++){
	// 	result_plaintexts[i] = new_TorusPolynomial(params->N);
	// }
	while(bound > 1e-5){
		std::cout << "Noise: " << bound << std::endl;
		for(int i = 0; i < t; i++){
			cycle_counts_partial[i] = 0;
			cycle_counts_final[i] = 0;
			for(int j = 0; j < params->N; j++){
				partial_ciphertexts[i]->coefsT[j] = 0;
			}
		}
		/* Each party performs partial decryption on its own */
		for(int i = 0; i < t; i++){
			partialDecrypt(ciphertext, params, partial_ciphertexts[i], cycle_counts_partial, i, subset, t, p, bound);
		}
		/* Each party performs final decryption on its own */
		for(int i = 0; i < t; i++){
			finalDecrypt(ciphertext, partial_ciphertexts, params, cycle_counts_final, i, subset, t, p);
		}
		bound /= 2;
	}
}