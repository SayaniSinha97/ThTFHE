# SubmissionCodes
Contains artifacts related to the submission

# Installation

```bash
git clone --recursive <repo url>
cd tfhe
make
sudo make install
rm -rf build
cd ..
cd OpenBLAS
make
sudo make PREFIX=/usr/local install
cd ..
```

Also install the Boost development library.

# User Installation

To install the asymmetric key version of the library, follow the Dev Installation setup
and then in the main directory perform `sudo make install`.

This installs `libthfhe` in `/usr/local/lib`.
To use this library with your own code, include `thfhe.hpp` in source code
and (statically) link the library with `-lthfhe` during compilation.

# Running code

Change `LD_LIBRARY_PATH`:

```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
```
Make:

```bash
make
```
bin/keygen generates keys.

```bash
./bin/keygen
```
bin/tlwetn encrypts an interger first to a Torus ring-LWE ciphertext, then does threshold secret sharing and then performs threshold decryption. To check (4,7)-threshold decryption by party 1, 3, 4, 5, run:

```bash
./bin/tlwetn 4 7 1 3 4 5
```
bin/convert takes two integers from test/plain22.txt and test/plain23.txt, encrypts them to Torus LWE ciphertexts, performs homomorphic AND, converts the resulting Torus LWE ciphertext into a Torus ring-LWE ciphertext and then performs threshold decryption on this ciphertext.

```bash
./bin/convert
```

bin/KNN_medical_data uses data from test/bootstrap_modules folder and implements our use-case.

```bash
./bin/KNN_medical_data
```

If a binary uses OpenMP and/or OpenBLAS, run it as:

```bash
OMP_NUM_THREADS=<cores> OPENBLAS_NUM_THREADS=<cores> ./<binary> <params>
```
The changes to `LD_LIBRARY_PATH`, `OMP_NUM_THREADS` and `OPENBLAS_NUM_THREADS` are to be done for using `libthfhe` in user code as well.

# Libthfhe sample program

The below program serves as an intended way to use `libthfhe` in user code.
It generates a public key and creates a (3, 5) threshold share of the secret key.
The public key is then used to encrypt two bits and there XOR is computed homomorphically.
Finally the threshold decryption of that ciphertext is done after converting it to TLWE.


```cpp
#include <tfhe/tfhe.h>
#include <thfhe.hpp>    // Main header needed

#include <iostream>
#include <string>
#include <fstream>
#include <cstdio>


int main()
{
    int t = 3;
    int T = 5;
    auto genKey = new ThFHE();
    genKey->KeyGen(t, T);
    auto params = initialize_gate_bootstrapping_params();

    LweSample *one = new_LweSample(params->in_out_params);
    LweSample *zero = new_LweSample(params->in_out_params);
    genKey->pk->Encrypt(one, 1);
    genKey->pk->Encrypt(zero, 0);

    LweSample *res = new_LweSample(params->in_out_params);
    bootsXOR(res, one, zero, genKey->bk);
    
    // Sanity Check by plain decryption
    std::cout << bootsSymDecrypt(res, genKey->sk) << std::endl;
    std::cout << bootsSymDecrypt(one, genKey->sk) << std::endl;
    std::cout << bootsSymDecrypt(zero, genKey->sk) << std::endl;

    ThFHEKeyShare *share1 = new ThFHEKeyShare();
    ThFHEKeyShare *share2 = new ThFHEKeyShare();
    ThFHEKeyShare *share3 = new ThFHEKeyShare();

    genKey->GetShareSet(1, share1);
    genKey->GetShareSet(2, share2);
    genKey->GetShareSet(3, share3);
    // In production, these shares are to be distributed to the parties
    // and the genKey object should be destroyed.

    TLweParams *tparams = new_TLweParams(1024, 1, pow(2., -25), pow(2, -15));

    TLweSample *tres = new_TLweSample(tparams);
    
    // This conversion is needed for the scheme to work
    // This has potential for ciphertext packing
    // Check the paper for details.
    TLweFromLwe(tres, res, tparams);
    
    TorusPolynomial **poly = new TorusPolynomial*[3];
    poly[0] = new_TorusPolynomial(1024);
    share1->PartialDecrypt(tres, tparams, poly[0], {1, 2, 3}, 3, 5, 0.0001);
    poly[1] = new_TorusPolynomial(1024);
    share2->PartialDecrypt(tres, tparams, poly[1], {1, 2, 3}, 3, 5, 0.0001);
    poly[2] = new_TorusPolynomial(1024);
    share3->PartialDecrypt(tres, tparams, poly[2], {1, 2, 3}, 3, 5, 0.0001);

    int msg = finalDecrypt(tres, poly, tparams, {1, 2, 3}, 3, 5);
    
    // This result must match with the sanity check results printed above
    std::cout << msg << std::endl;


    return 0;
}
```


# Running on Raspberry Pi

Follow the README.md file in `rpi` directory.
