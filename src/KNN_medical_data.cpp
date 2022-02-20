#include <iostream>
#include<bits/stdc++.h>
#include<fstream>
#include <sstream>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <time.h>

#define MAXLEN 4

#include <tfhe/lwe-functions.h>
#include <tfhe/numeric_functions.h>
#include <tfhe/tlwe_functions.h>
#include <random>
#include "threshold_decryption_functions.hpp"







void keygen(){

    const int minimum_lambda = 110;
    auto params = new_default_gate_bootstrapping_parameters(minimum_lambda);

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


}



void encrypt(){

    FILE *skFile = fopen("test/secret.key", "rb");
    auto key_ = new_tfheGateBootstrappingSecretKeySet_fromFile(skFile);
    fclose(skFile);

    FILE *paramFile = fopen("test/secret.params", "rb");
    auto params_ = new_tfheGateBootstrappingParameterSet_fromFile(paramFile);
    fclose(paramFile); 

    LweSample *allOne = new_gate_bootstrapping_ciphertext_array(32, params_);
    LweSample *allZero = new_gate_bootstrapping_ciphertext_array(32, params_);
    LweSample *lsbOne = new_gate_bootstrapping_ciphertext_array(32, params_);
    LweSample *lsbZero = new_gate_bootstrapping_ciphertext_array(32, params_);

    for (int i = 0; i < 32; i++){
        bootsSymEncrypt(&allOne[31-i], 1, key_);
    }

    for (int i = 0; i < 32; i++){
        bootsSymEncrypt(&allZero[31-i], 0, key_);
    }

    for (int i = 0; i < 32; i++){
        if(i==0)
            bootsSymEncrypt(&lsbOne[31-i], 1, key_);
        else 
            bootsSymEncrypt(&lsbOne[31-i], 0, key_);
    }

    for (int i = 0; i < 32; i++){
        if(i==0)
            bootsSymEncrypt(&lsbZero[31-i], 0, key_);
        else 
            bootsSymEncrypt(&lsbZero[31-i], 1, key_);
    }

    FILE *all_one = fopen("test/bootstrap_modules/allOne.data", "wb");
    for (int i = 0; i < 32; i++){
        export_gate_bootstrapping_ciphertext_toFile(all_one, &allOne[i], params_);
    }
    fclose(all_one);
    delete_gate_bootstrapping_ciphertext_array(32, allOne);

    FILE *all_zero = fopen("test/bootstrap_modules/allZero.data", "wb");
    for (int i = 0; i < 32; i++){
        export_gate_bootstrapping_ciphertext_toFile(all_zero, &allZero[i], params_);
    }
    fclose(all_zero);
    delete_gate_bootstrapping_ciphertext_array(32, allZero);


    FILE *lsb_one = fopen("test/bootstrap_modules/lsbOne.data", "wb");
    for (int i = 0; i < 32; i++){
        export_gate_bootstrapping_ciphertext_toFile(lsb_one, &lsbOne[i], params_);
    }
    fclose(lsb_one);
    delete_gate_bootstrapping_ciphertext_array(32, lsbOne);

    FILE *lsb_zero = fopen("test/bootstrap_modules/lsbZero.data", "wb");
    for (int i = 0; i < 32; i++){
        export_gate_bootstrapping_ciphertext_toFile(lsb_zero, &lsbZero[i], params_);
    }
    fclose(lsb_zero);
    delete_gate_bootstrapping_ciphertext_array(32, lsbZero);


    delete_gate_bootstrapping_parameters(params_);

}





void onesComp(LweSample *result, LweSample *input1, LweSample *input2, const int nbits, const TFheGateBootstrappingCloudKeySet *bk)
{
    for (int i = 0; i < nbits; i++){
        bootsXOR(&result[i], &input1[i], &input2[i], bk);
    }
}

void FullAdder(LweSample *sum2, LweSample *carrybit, LweSample *input1, LweSample *input2, const int nbits, const TFheGateBootstrappingCloudKeySet *bk)
{
    LweSample *sum1 = new_gate_bootstrapping_ciphertext_array(32, bk->params);
    LweSample *carry1 = new_gate_bootstrapping_ciphertext_array(32, bk->params);
    LweSample *carry2 = new_gate_bootstrapping_ciphertext_array(32, bk->params);

    for (int i = nbits-1; i >= 0; i--){
        // half adder 1
        bootsXOR(&sum1[i], &input1[i], &input2[i], bk);
        bootsAND(&carry1[i], &input1[i], &input2[i], bk);

        // half adder 2
        bootsXOR(&sum2[i], &sum1[i], &carrybit[i], bk);
        bootsAND(&carry2[i], &sum1[i], &carrybit[i], bk);

        // final carry
        if(i != 0)
        bootsOR(&carrybit[i-1], &carry1[i], &carry2[i], bk);
    }

    delete_gate_bootstrapping_ciphertext_array(32, sum1);
    delete_gate_bootstrapping_ciphertext_array(32, carry1);
    delete_gate_bootstrapping_ciphertext_array(32, carry2);
}



void difference(LweSample *diff, LweSample *ciphertext_1, LweSample *ciphertext_2, const int nbits ){
    FILE *cloudKeyFile = fopen("test/cloud.key", "rb");
    auto cloud_key_ = new_tfheGateBootstrappingCloudKeySet_fromFile(cloudKeyFile);
    fclose(cloudKeyFile);

    // 2's complement of ciphertext 2
    LweSample *onesComplement = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
    LweSample *twosComplement = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);

    LweSample *allOne = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
    FILE *all_one = fopen("test/bootstrap_modules/allOne.data", "rb");
    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(all_one, &allOne[i], cloud_key_->params);
    fclose(all_one);

    LweSample *lsbOne = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
    FILE *lsb_one = fopen("test/bootstrap_modules/lsbOne.data", "rb");
    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(lsb_one, &lsbOne[i], cloud_key_->params);
    fclose(lsb_one);

    LweSample *carry2scomplement = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
    FILE *lsb_zero = fopen("test/bootstrap_modules/lsbZero.data", "rb");
    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(lsb_zero, &carry2scomplement[i], cloud_key_->params);
    fclose(lsb_zero);

    onesComp(onesComplement, allOne, ciphertext_2, 32, cloud_key_);

    // Two's complement of ciphertext2
    FullAdder(twosComplement, carry2scomplement, onesComplement, lsbOne, 32, cloud_key_);


    //LweSample *difference = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
    LweSample *borrow = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);

    FILE *f = fopen("test/bootstrap_modules/lsbZero.data", "rb");
    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(f, &borrow[i], cloud_key_->params);
    fclose(f);

    //diff = cipheretext1 - ciphertext2
    FullAdder(diff, borrow, ciphertext_1, twosComplement, 32, cloud_key_);


    delete_gate_bootstrapping_ciphertext_array(32, onesComplement);
    delete_gate_bootstrapping_ciphertext_array(32, twosComplement);
    delete_gate_bootstrapping_ciphertext_array(32, allOne);
    delete_gate_bootstrapping_ciphertext_array(32, lsbOne);
    delete_gate_bootstrapping_ciphertext_array(32, borrow);
    delete_gate_bootstrapping_ciphertext_array(32, carry2scomplement);
    delete_gate_bootstrapping_cloud_keyset(cloud_key_);
}



void distance(LweSample *dist, LweSample *input1, LweSample *input2, const int n, const TFheGateBootstrappingCloudKeySet *cloud_key_){

    LweSample *difference1 = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
    LweSample *difference2 = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);

    difference(difference1, input1, input2, n);
    difference(difference2, input2, input1, n);

    //Homomorphic bootstrapped Mux(a,b,c) = a?b:c = a*b + not(a)*c
    for(int i=0; i<n; i++){
        bootsMUX(&dist[i], &difference1[0], &difference2[i], &difference1[i], cloud_key_);
    }

    // printf("\n------------Calculating Distance ------\n");
    // int dis =  decrypt_cipher(dist);
    // printf("Distance = %d\n", dis);

    delete_gate_bootstrapping_ciphertext_array(MAXLEN*8, difference1);
    delete_gate_bootstrapping_ciphertext_array(MAXLEN*8, difference2);
}


void distance_bw_data(LweSample *result, LweSample **input1, LweSample **input2, int col_size, int nbits,  const TFheGateBootstrappingCloudKeySet *cloud_key_){
    LweSample *dist = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);

    //LweSample *man_dist = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
    FILE *all_zero = fopen("test/bootstrap_modules/allZero.data", "rb");
    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(all_zero, &result[i], cloud_key_->params);
    fclose(all_zero);

    LweSample *carry2 = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
    FILE *lsb_zero = fopen("test/bootstrap_modules/lsbZero.data", "rb");
    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(lsb_zero, &carry2[i], cloud_key_->params);
    fclose(lsb_zero);

    for(int i=1; i<col_size; i++){
        distance(dist, input1[i], input2[i], nbits, cloud_key_);
        FullAdder(result, carry2, result, dist, nbits, cloud_key_);
    }

    printf("\n------------Computing Distance ------\n");
    //int dis =  decrypt_cipher(result);
    //printf("Total Distance = %d\n", dis);
    delete_gate_bootstrapping_ciphertext_array(32, carry2);
}


void encrypt_data(LweSample *cipher, int data){
    FILE *skFile = fopen("test/secret.key", "rb");
    auto key_ = new_tfheGateBootstrappingSecretKeySet_fromFile(skFile);
    fclose(skFile);

    FILE *paramFile = fopen("test/secret.params", "rb");
    auto params_ = new_tfheGateBootstrappingParameterSet_fromFile(paramFile);
    fclose(paramFile); 

    for(int j = 0; j < 32; j++){
        bootsSymEncrypt(&cipher[31-j], (data >> j) & 1, key_);
    }

    delete_gate_bootstrapping_secret_keyset(key_);
    delete_gate_bootstrapping_parameters(params_);
}

void encrypt_dataset(LweSample **ciphers, std::vector<int>&row){

    FILE *skFile = fopen("test/secret.key", "rb");
    auto key_ = new_tfheGateBootstrappingSecretKeySet_fromFile(skFile);
    fclose(skFile);

    FILE *paramFile = fopen("test/secret.params", "rb");
    auto params_ = new_tfheGateBootstrappingParameterSet_fromFile(paramFile);
    fclose(paramFile); 

    for(int i=0; i<row.size();i++){

        for(int j = 0; j < 32; j++){
            bootsSymEncrypt(&ciphers[i][31-j], (row[i] >> j) & 1, key_);
        }
    }

    delete_gate_bootstrapping_secret_keyset(key_);
    delete_gate_bootstrapping_parameters(params_);
}

int decrypt_bit(LweSample *cipher){
    FILE *skFile = fopen("test/secret.key", "rb");
    auto key_ = new_tfheGateBootstrappingSecretKeySet_fromFile(skFile);
    fclose(skFile);

    printf("........Decrypting cipher......\n");

    int num = 0;
    int bit = bootsSymDecrypt(&cipher[0], key_);
    num = num*2 + bit;
    printf(" %d\n",num);
    printf("\n Decryption over --\n");
    return num;
    
}


int decrypt_cipher(LweSample *cipher){
    FILE *skFile = fopen("test/secret.key", "rb");
    auto key_ = new_tfheGateBootstrappingSecretKeySet_fromFile(skFile);
    fclose(skFile);

    printf("........Decrypting cipher......\n");

    int num = 0;
    for (int j = 0; j < 32; j++){
        int bit = bootsSymDecrypt(&cipher[j], key_);
        num = num*2 + bit;
    }
    printf(" %d\n",num);
    printf("\n Decryption over --\n");
    return num;
    
}


void decrypt_ciphers(LweSample **cipher, int n){
    FILE *skFile = fopen("test/secret.key", "rb");
    auto key_ = new_tfheGateBootstrappingSecretKeySet_fromFile(skFile);
    fclose(skFile);

    printf("........Decrypting ciphers......\n");

    for(int i=0;i<n;i++){
        int num = 0;
        for (int j = 0; j < 32; j++){
            int bit = bootsSymDecrypt(&cipher[i][j], key_);
            num = num*2 + bit;
        }
        printf(" %d", num);
    }
    printf("\n Decryption Done \n");

    delete_gate_bootstrapping_secret_keyset(key_);
}



void bubble_sort(LweSample **cipher, int n){

    printf("\n Inside Bubble Sort \n");

    FILE *cloudKeyFile = fopen("test/cloud.key", "rb");
    auto cloud_key_ = new_tfheGateBootstrappingCloudKeySet_fromFile(cloudKeyFile);
    fclose(cloudKeyFile);

    LweSample *allZero = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
    FILE *all_zero = fopen("test/bootstrap_modules/allZero.data", "rb");
    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(all_zero, &allZero[i], cloud_key_->params);
    fclose(all_zero);


    LweSample *diff = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
    LweSample *ans1= new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
    LweSample *ans2 = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);

    int count = n;
    while(count--){
        for(int i=1; i<n; i++){
            difference(diff, cipher[i-1], cipher[i], 32);

            for(int j=0;j<32;j++){
            // bigger element
                bootsMUX(&ans1[j],  &diff[0], &cipher[i][j], &cipher[i-1][j], cloud_key_);
            
            // smaller element
                bootsMUX(&ans2[j],  &diff[0], &cipher[i-1][j], &cipher[i][j], cloud_key_);
            }

            for(int j=0; j<32;j++){
                bootsXOR(&cipher[i-1][j], &ans2[j], &allZero[j], cloud_key_);
                bootsXOR(&cipher[i][j], &ans1[j], &allZero[j], cloud_key_);
            }
        }
    }
    delete_gate_bootstrapping_ciphertext_array(32, allZero);
    delete_gate_bootstrapping_ciphertext_array(32, diff);
    delete_gate_bootstrapping_ciphertext_array(32, ans1);
    delete_gate_bootstrapping_ciphertext_array(32, ans2);
    delete_gate_bootstrapping_cloud_keyset(cloud_key_);

}



void sort_with_distance(LweSample ***sorted_train_data, LweSample **distant, int n, int col_size,  const TFheGateBootstrappingCloudKeySet *cloud_key_){
    LweSample *allZero = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
    FILE *all_zero = fopen("test/bootstrap_modules/allZero.data", "rb");
    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(all_zero, &allZero[i], cloud_key_->params);
    fclose(all_zero);


    LweSample *diff = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
    LweSample *ans1= new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
    LweSample *ans2 = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);

    LweSample **array1 = new LweSample*[col_size];
    LweSample **array2 = new LweSample*[col_size];

    for(int i=0;i<col_size;i++){
        array1[i] = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
        array2[i] = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
    }

    printf("\n ----- Inside sort by disatance1 ----- \n");
    int count = n;
    while(count--){

        for(int i=1; i<n; i++){
            difference(diff, distant[i-1], distant[i], 32);

            for(int j=0;j<32;j++){
                // bigger element
                bootsMUX(&ans1[j],  &diff[0], &distant[i][j], &distant[i-1][j], cloud_key_);

                // smaller element
                bootsMUX(&ans2[j],  &diff[0], &distant[i-1][j], &distant[i][j], cloud_key_);
            }
            // printf("\n Smaller = %d \n",decrypt_cipher(ans2));
            // printf("\n Bigger = %d \n",decrypt_cipher(ans1));

            for(int k=0; k<col_size; k++){
                for(int j =0; j<32; j++){
                    bootsMUX(&array1[k][j], &diff[0], &sorted_train_data[i][k][j], &sorted_train_data[i-1][k][j], cloud_key_);

                    bootsMUX(&array2[k][j], &diff[0], &sorted_train_data[i-1][k][j], &sorted_train_data[i][k][j], cloud_key_);
                }
            }

            // printf("\n ----- Inside sort by disatance3 ----- \n");
            // decrypt_ciphers(array2, col_size);
            // decrypt_ciphers(array1, col_size);
            

            for(int j=0; j<32;j++){
            //printf(" %d",j);
                bootsXOR(&distant[i-1][j], &ans2[j], &allZero[j], cloud_key_);
                bootsXOR(&distant[i][j], &ans1[j], &allZero[j], cloud_key_);

            }
            // printf("\n ----- Inside sort by disatance4 ----- \n");
            for(int k =0; k<col_size; k++){
                for(int j = 0; j<32; j++){
                //printf("k=%d j=%d\n", k, j);
                    bootsXOR(&sorted_train_data[i-1][k][j], &array2[k][j], &allZero[j], cloud_key_);
                    bootsXOR(&sorted_train_data[i][k][j], &array1[k][j], &allZero[j], cloud_key_);
                }
            }

        // Checking the claculation
        // for(int i=0; i<n;i++){
        //     decrypt_ciphers(sorted_train_data[i], col_size);
        //     decrypt_cipher(distant[i]);
        // }

        }

    }

        delete_gate_bootstrapping_ciphertext_array(32, allZero);
        delete_gate_bootstrapping_ciphertext_array(32, diff);
        delete_gate_bootstrapping_ciphertext_array(32, ans1);
        delete_gate_bootstrapping_ciphertext_array(32, ans2);
}



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


void ciphertext_conversion_threshold_decryption(LweSample *decisional_bit){
    FILE *skFile = fopen("test/secret.key", "rb");
    auto key = new_tfheGateBootstrappingSecretKeySet_fromFile(skFile);
    fclose(skFile);
    
    FILE *paramFile = fopen("test/secret.params", "rb");
    auto params = new_tfheGateBootstrappingParameterSet_fromFile(paramFile);
    fclose(paramFile);


    TLweParams *tlwe_params = new_TLweParams(1024, 1, 0.01, 0.2);
    TLweKey *tlwe_key = new_TLweKey(tlwe_params);
	tLweKeyGen(tlwe_key);

	TLweSample *resultOfEvalT;
	TorusPolynomial *result_plaintext;

	TLweKeyFromLweKey(key->lwe_key, tlwe_key);

    shareSecret(3, 5, tlwe_key, tlwe_params);

    printf("\n Threshold decryption on TLWE ciphertext  \n");

    double bound = 0.0125;
    int result_msg;
    std::vector<int> subset{1,2,4};
    while(bound > 1e-3){
	    result_msg = 0;
	    //Assuming k = 1, n = N = 1024, and converting lwe ciphertext of each of the result bit into corresponding ring-lwe ciphertext one by one and decrypting
	    //TODO: Pack all 32 lwe ciphertexts into one tlwe ciphertext and call threshold_decryption once
	    //NEXT TODO: Do the same but with k > 1

		resultOfEvalT = new_TLweSample(tlwe_params);
		TLweFromLwe(resultOfEvalT, &decisional_bit[0], tlwe_params);
		result_plaintext = new_TorusPolynomial(tlwe_params->N);
		thresholdDecrypt(result_plaintext, resultOfEvalT, tlwe_params, subset, 3, 5, bound);
		result_msg += (result_plaintext->coefsT[0] > 0 ? 1 : 0);

		std::cout << "bound: " << bound << "   result_msg: " << result_msg << "\n";
		bound /= 2;
	}
}



void inputDataSet()
{
    FILE *cloudKeyFile = fopen("test/cloud.key", "rb");
    auto cloud_key_ = new_tfheGateBootstrappingCloudKeySet_fromFile(cloudKeyFile);
    fclose(cloudKeyFile);

    int col_size = 14;
    int train_row_size = 5;
    int test_row_size = 1;
    std::vector<std::vector<int>> row(train_row_size + test_row_size);
    std::string line, word, temp;

    std::ifstream read("test/bootstrap_modules/data1.csv"); 
    read>>line;

    LweSample ***cipher_train_dataset = new LweSample**[train_row_size];
    LweSample ***cipher_test_dataset = new LweSample**[test_row_size];

    for(int i=0;i<train_row_size;i++){
        cipher_train_dataset[i] = new LweSample*[col_size];
        for(int j=0; j<col_size; j++){
            cipher_train_dataset[i][j] = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
        }   
    }

    for(int i=0;i<test_row_size;i++){
        cipher_test_dataset[i] = new LweSample*[col_size];
        for(int j=0; j<col_size; j++){
            cipher_test_dataset[i][j] = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
        }   
    }

    int row_count = 0;
    while (read >> line)
    {
        row[row_count] = std::vector<int>(col_size);
        row[row_count].clear();
        //std::cout << line << std::endl;

        std::stringstream s(line);
  
        while (std::getline(s, word, ',')) {
        //std::cout << word<<" ";

        std::stringstream ss(word);
        int x =0;
        ss >> x;
        row[row_count].push_back(x);
        }
        //std::cout<<"\n";

        if(row_count < train_row_size){
            encrypt_dataset(cipher_train_dataset[row_count], row[row_count]);
            decrypt_ciphers(cipher_train_dataset[row_count], col_size);
        }
        else{
            encrypt_dataset(cipher_test_dataset[row_count-train_row_size], row[row_count]);
            decrypt_ciphers(cipher_test_dataset[row_count-train_row_size], col_size);
        }

        //encrypt_dataset(cipher_dataset[row_count], row[row_count]);
        //distance(cipher_dataset[row_count][1], cipher_dataset[row_count][2], 32);
        //decrypt_ciphers(cipher_dataset[row_count], col_size);
        //bubble_sort(cipher_dataset[row_count], col_size);
        //decrypt_ciphers(cipher_dataset[row_count], col_size);


    row_count++;
    printf(" Row #%d\n",row_count);
    if(row_count == (train_row_size + test_row_size))
        break;
    }


    int K = 5;
    int correct = 0;
    int incorrect = 0;

    int threshold = K/2;
    LweSample *threshold_cipher = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
    encrypt_data(threshold_cipher, threshold);
    //decrypt_cipher(threshold_cipher);

    LweSample *allOne = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
    FILE *all_one = fopen("test/bootstrap_modules/allOne.data", "rb");
    for (int i = 0; i < 32; i++)
        import_gate_bootstrapping_ciphertext_fromFile(all_one, &allOne[i], cloud_key_->params);
    fclose(all_one);


    LweSample **distant = new LweSample*[train_row_size];
    LweSample ***sorted_train_dataset = new LweSample**[train_row_size];



    for(int i=0;i<train_row_size;i++){
        distant[i] = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
        sorted_train_dataset[i] = new LweSample*[col_size];
        for(int j=0; j<col_size; j++){
            sorted_train_dataset[i][j] = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
        }
    }
//    ------ Testing Started -----
    for(int i=0; i<test_row_size;i++){
        printf("\n------ Test data #%d -------- \n",i+1);
        #pragma omp parallel for
            for(int j=0; j<train_row_size; j++){
                distance_bw_data(distant[j], cipher_test_dataset[i], cipher_train_dataset[j], col_size-1, 32, cloud_key_);

            // Copying in sorted_train_dataset
                for(int k=0; k<col_size; k++){
                    for(int l=0; l<32;l++){
                        bootsMUX(&sorted_train_dataset[j][k][l], &allOne[l],  &cipher_train_dataset[j][k][l], &cipher_train_dataset[j][k][l], cloud_key_ );
                    }
                }
            }
        for(int j=0; j<train_row_size;j++){
            decrypt_cipher(distant[j]);
        }

        //bubble_sort(distant, train_row_size);

        // Sorting train data  distance-wise
        sort_with_distance(sorted_train_dataset, distant, train_row_size, col_size, cloud_key_);

        // printf("\n ---- Checking sorted train data ----\n");
        // for(int j=0; j<train_row_size;j++){
        //     decrypt_ciphers(sorted_train_dataset[j], col_size);
        //     decrypt_cipher(distant[j]);
        // }
        LweSample *count = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
        FILE *all_zero = fopen("test/bootstrap_modules/allZero.data", "rb");
        for (int i = 0; i < 32; i++)
            import_gate_bootstrapping_ciphertext_fromFile(all_zero, &count[i], cloud_key_->params);
        fclose(all_zero);

        LweSample *carry2 = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
        FILE *lsb_zero = fopen("test/bootstrap_modules/lsbZero.data", "rb");
        for (int i = 0; i < 32; i++)
            import_gate_bootstrapping_ciphertext_fromFile(lsb_zero, &carry2[i], cloud_key_->params);
        fclose(lsb_zero);

        printf("\n Count reset to  %d\n", decrypt_cipher(count));
        for(int j = 0; j < K; j++){
            FullAdder(count, carry2, count, sorted_train_dataset[j][col_size-1], 32, cloud_key_ );
        }

        LweSample *decisional_bit = new_gate_bootstrapping_ciphertext_array(1, cloud_key_->params);
        LweSample *diff = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
        difference(diff, threshold_cipher, count, 32);

        LweSample *allZero = new_gate_bootstrapping_ciphertext_array(32, cloud_key_->params);
        all_zero = fopen("test/bootstrap_modules/allZero.data", "rb");
        for (int i = 0; i < 32; i++)
            import_gate_bootstrapping_ciphertext_fromFile(all_zero, &allZero[i], cloud_key_->params);
        fclose(all_zero);
        bootsXOR(&decisional_bit[0], &diff[0], &allZero[0], cloud_key_ );

        // int p_count = decrypt_cipher(count);
        int original_decision = decrypt_cipher(cipher_test_dataset[i][col_size-1]);

        // given as input to harware decryption algorithm
        int predicted_bit = decrypt_bit(decisional_bit);

        if(predicted_bit == original_decision){
            printf("\n Correctly predicted \n");
            correct++; 
        }
        else{
            printf("\n incorrect prediction \n");
            incorrect++;    
        }


    ciphertext_conversion_threshold_decryption(decisional_bit);

    delete_gate_bootstrapping_ciphertext_array(MAXLEN*8, count);
    delete_gate_bootstrapping_ciphertext_array(MAXLEN*8, carry2);
    delete_gate_bootstrapping_ciphertext_array(MAXLEN*8, allZero);
    delete_gate_bootstrapping_ciphertext_array(1, decisional_bit);
    }

    float acc = (float)correct/(correct+incorrect);
    printf("\n Accuracy = %f \n", acc);


    delete_gate_bootstrapping_cloud_keyset(cloud_key_);

    for(int i=0; i<train_row_size; i++){
        for(int j=0; j<col_size; j++){
            delete_gate_bootstrapping_ciphertext_array(MAXLEN*8, cipher_train_dataset[i][j]);
        }
    }

    for(int i=0; i<test_row_size; i++){
        for(int j=0; j<col_size; j++){
            delete_gate_bootstrapping_ciphertext_array(MAXLEN*8, cipher_test_dataset[i][j]);
        }
    }
        
}




void test_lwe2rlwe_conversion(){
    FILE *skFile = fopen("test/secret.key", "rb");
    auto key_ = new_tfheGateBootstrappingSecretKeySet_fromFile(skFile);
    fclose(skFile);

    FILE *paramFile = fopen("test/secret.params", "rb");
    auto params_ = new_tfheGateBootstrappingParameterSet_fromFile(paramFile);
    fclose(paramFile);    

    int len;
    int buff1=1;

    LweSample *ciphertext1 = new_gate_bootstrapping_ciphertext_array(1, params_);
    printf("\n------Inside Encrypt bit function----\n");

    bootsSymEncrypt(&ciphertext1[0], buff1, key_);
    // int bit = bootsSymDecrypt(&ciphertext1[0], key_);
    // printf("%d ", bit);

    ciphertext_conversion_threshold_decryption(ciphertext1);


    delete_gate_bootstrapping_ciphertext_array(1, ciphertext1);
    delete_gate_bootstrapping_secret_keyset(key_);
    delete_gate_bootstrapping_parameters(params_);
}



float big_job(int x){
    return x+x;
}

void consume(float x, float res){
    printf("Consume called \n");
}

int main()
{
    clock_t start1, end1;
    start1 = clock();

    keygen();
    encrypt();

    // float res;
    // #pragma omp parallel for
    // { 
    //     float B; int i, id, nthrds;
    //     id = omp_get_thread_num();
    //     nthrds = omp_get_num_threads();

    //     for(i=id;i<20;i++){
    //     B = big_job(i);
    //     printf("Value = %f \n",B);
    //     #pragma omp critical
    //     consume (B, res);
    //     }
    // }


    //test_lwe2rlwe_conversion();

    inputDataSet();

    end1 = clock();

    printf("\nTotal time : %.3f ms \n", (double(end1-start1)/CLOCKS_PER_SEC)*1000);

return 0;
}