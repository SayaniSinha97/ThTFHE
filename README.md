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

Modify Compute function in src/Compute.cpp to create new circuits.
Modify test.py accordingly.

bin/keygen generates keys.
bin/encrypt takes data from test/plain.txt and encrypts it.
bin/compute computes on the encryptions.
bin/decrypt performs the decryption.

If a binary uses OpenMP and/or OpenBLAS, run it as:

```bash
OMP_NUM_THREADS=<cores> OPENBLAS_NUM_THREADS=<cores> ./<binary> <params>
```

# Running on Raspberry Pi

Follow the README.md file in `rpi` directory.
