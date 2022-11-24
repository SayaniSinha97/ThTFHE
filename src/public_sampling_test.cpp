#include <iostream>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <tfhe/lwe-functions.h>
#include <tfhe/numeric_functions.h>
#include <tfhe/tfhe_garbage_collector.h>
#include <random>
#include <time.h>

typedef struct myParameterSet{
	TLweParams *tlweparams;
	TFheGateBootstrappingParameterSet *bootparams;
}myParameterSet;
myParameterSet *initialize_gate_bootstrapping_params() {
    static const int32_t N = 1024;
    static const int32_t k = 1;
    static const int32_t n = 1024;
    static const int32_t bk_l = 3;
    static const int32_t bk_Bgbit = 7;
    static const int32_t ks_basebit = 2;
    static const int32_t ks_length = 8;
    static const double ks_stdev = pow(2.,-15); //standard deviation
    static const double bk_stdev = pow(2.,-25);; //standard deviation
    static const double max_stdev = 0.012467; //max standard deviation for a 1/4 msg space

    LweParams *params_in = new_LweParams(n, ks_stdev, max_stdev);
    TLweParams *params_accum = new_TLweParams(N, k, bk_stdev, max_stdev);
    TGswParams *params_bk = new_TGswParams(bk_l, bk_Bgbit, params_accum);

    TfheGarbageCollector::register_param(params_in);
    TfheGarbageCollector::register_param(params_accum);
    TfheGarbageCollector::register_param(params_bk);

    auto bootparams = new TFheGateBootstrappingParameterSet(ks_length, ks_basebit, params_in, params_bk);
    myParameterSet *params = new myParameterSet();
    params->tlweparams = params_accum;
    params->bootparams = bootparams;
    return params;
}

void BitwiseEncrypt(LweSample *result_array, int32_t msg, TFheGateBootstrappingSecretKeySet *key){
	for (int i = 0; i < 32; i++){
        bootsSymEncrypt(&result_array[i], (msg >> i) & 1, key);
    }
}

int BitwiseDecrypt(LweSample *cipher, TFheGateBootstrappingSecretKeySet *key){
	int res = 0;
	for(int i = 0; i < 32; i++){
		res += bootsSymDecrypt(&cipher[i], key) << i;
	}
	return res;
}

void add_msg(LweSample *result, const LweParams *lweparams, int32_t msg){
	int n = lweparams->n;
	for(int i = 0; i < 32; i++){
		(&result[i])->b += modSwitchToTorus32(((msg >> i) & 1), 2);
	}
}

void generate_sample(LweSample *sample, LweSample *input, const LweParams *lweparams, const TFheGateBootstrappingCloudKeySet *cloud_key, int bit, int32_t msg0, int32_t msg1){
	for(int i = 0; i < 32; i++){
		bootsXOR(&sample[i], &input[i], &input[i], cloud_key);
	}

	if(bit == 0){
		add_msg(sample, lweparams, msg0);
	}
	else{
		add_msg(sample, lweparams, msg1);
	}
}

int main(){
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dist(0, 1);

	auto myparams = initialize_gate_bootstrapping_params();
	auto params = myparams->bootparams;
	const LweParams *lweparams = params->in_out_params;

	TFheGateBootstrappingSecretKeySet *key = new_random_gate_bootstrapping_secret_keyset(params);
	auto cloud_key = &key->cloud;

	/* Read two different plaintexts msg0 and msg1 */
	FILE *plaintext0 = fopen("test/m0.txt", "r");
    int32_t msg0;
    fscanf(plaintext0, "%d", &msg0);
    fclose(plaintext0);
    std::cout << "Plaintext0: " << msg0 << "\n";

    FILE *plaintext1 = fopen("test/m1.txt", "r");
    int32_t msg1;
    fscanf(plaintext1, "%d", &msg1);
    fclose(plaintext1);
    std::cout << "Plaintext1: " << msg1 << "\n";

	/* choose a random bit */
	int sel = dist(gen);

	/* Encryption : Encrypt message m0 if sel = 0, else encrypt message m1 */
	LweSample *ciphertext = new_gate_bootstrapping_ciphertext_array(32, params);
    
	if(sel == 0){		
		BitwiseEncrypt(ciphertext, msg0, key);
	}
	else{
		BitwiseEncrypt(ciphertext, msg1, key);
	}
	std::cout << "passing encryption of plaintext: " << (sel == 0 ? msg0 : msg1) << "\n\n";

	LweSample *fresh_sample0 = new_gate_bootstrapping_ciphertext_array(32, params);
	LweSample *fresh_sample1 = new_gate_bootstrapping_ciphertext_array(32, params);

	/* Generate fresh samples from two distributions and verify correctness */
	std::cout << "Generating sample of encryption of " << msg0 << "\n";
	generate_sample(fresh_sample0, ciphertext, lweparams, cloud_key, 0, msg0, msg1);
	int32_t dmsg0 = BitwiseDecrypt(fresh_sample0, key);
	std::cout << "expected: " << msg0 << ", generated: " << dmsg0 << "\n\n";

	std::cout << "Generating sample of encryption of " << msg1 << "\n";
	generate_sample(fresh_sample1, ciphertext, lweparams, cloud_key, 1, msg0, msg1);
	int32_t dmsg1 = BitwiseDecrypt(fresh_sample1, key);
	std::cout << "expected: " << msg1 << ", generated: " << dmsg1 << "\n\n";
	
	return 0;
}