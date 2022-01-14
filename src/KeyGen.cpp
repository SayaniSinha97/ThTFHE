#include <iostream>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <tfhe/tfhe_garbage_collector.h>


TFheGateBootstrappingParameterSet *initialize_gate_bootstrapping_params() {
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

    return new TFheGateBootstrappingParameterSet(ks_length, ks_basebit, params_in, params_bk);
}


int main()
{
    // const int minimum_lambda = 110;
    auto params = initialize_gate_bootstrapping_params();

    uint32_t seed[] = { 100, 20032, 21341 };
    tfhe_random_generator_setSeed(seed, 3);

    auto key = new_random_gate_bootstrapping_secret_keyset(params);

    FILE* secret_key = fopen("test/secret.key", "wb");
    export_tfheGateBootstrappingSecretKeySet_toFile(secret_key, key);
    fclose(secret_key);

    FILE* cloud_key = fopen("test/cloud.key", "wb");
    export_tfheGateBootstrappingCloudKeySet_toFile(cloud_key, &key->cloud);
    fclose(cloud_key);

    FILE* secret_params = fopen("test/secret.params", "wb");
    export_tfheGateBootstrappingParameterSet_toFile(secret_params, params);
    fclose(secret_params);

    delete_gate_bootstrapping_secret_keyset(key);
    delete_gate_bootstrapping_parameters(params);

    return 0;
}
