# Running on Raspbian Lite

Make sure to install the packages `raspberrypi-kernel` and `raspberrypi-kernel-headers`.
Then do the following:

```bash
cd src/utils
make
sudo insmod enable_ccr.ko
```

This will enable cycle counter register reading from userspace.

Then run `make -e TP={ 47 | 34 | 35 | 59 | 610 | 711 | 913 | 1015 }`
to build for 4 out of 7, 3 out of 4 and so on threshold decryption systems.
Then run the `main` binary.

This package has the test data included in `include/tfhe/rlwe_test_data*.h` files.