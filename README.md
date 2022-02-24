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
bin/convert takes two integers from test/plain22.txt and test/plain23.txt, encrypts them to Torus LWE ciphertexts, performs homomorphic AND, convert the resulting Torus LWE ciphertext into a Torus ring-LWE ciphertext and then perform threshold decryption on this ciphertext.

```bash
./bin/convert
```

If a binary uses OpenMP and/or OpenBLAS, run it as:

```bash
OMP_NUM_THREADS=<cores> OPENBLAS_NUM_THREADS=<cores> ./<binary> <params>
```

# Running on Raspberry Pi

Follow the README.md file in `rpi` directory.
