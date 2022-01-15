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

bin/keygen generates keys.
bin/tlwetn encrypts an interger first, then does threshold secret sharing and then performs threshold decryption.
bin/convert takes two integers, encrypts them, performs homomorphic AND, and then threshold decrypt the resulting ciphertext.

If a binary uses OpenMP and/or OpenBLAS, run it as:

```bash
OMP_NUM_THREADS=<cores> OPENBLAS_NUM_THREADS=<cores> ./<binary> <params>
```

# Running on Raspberry Pi

Follow the README.md file in `rpi` directory.
