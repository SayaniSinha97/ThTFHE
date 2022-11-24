TFHE_PREFIX = /usr/local
TFHE_OPTIONS = C_INCLUDE_PATH=$$C_INCLUDE_PATH:${TFHE_PREFIX}/include CPLUS_INCLUDE_PATH=$$CPLUS_INCLUDE_PATH:${TFHE_PREFIX}/include LIBRARY_PATH=$$LIBRARY_PATH:${TFHE_PREFIX}/lib LD_LIBRARY_PATH=$$LD_LIBRARY_PATH:${TFHE_PREFIX}/lib
TFHE_LIB = -ltfhe-spqlios-fma
BLAS_LIB = -lopenblas

all: bin/keygen bin/tlwetn bin/convert bin/KNN_medical_data bin/libthfhe.a bin/public_sampling_test

bin/public_sampling_test: src/public_sampling_test.cpp
	${TFHE_OPTIONS} g++ -g src/public_sampling_test.cpp -o bin/public_sampling_test ${TFHE_LIB}

bin/libthfhe.a: src/libthfhe.cpp src/thfhe.hpp
	rm -f bin/libthfhe.a
	${TFHE_OPTIONS} g++ -g -fopenmp -O3 -c src/libthfhe.cpp -o bin/libthfhe.o ${TFHE_LIB} ${BLAS_LIB}
	ar -crs bin/libthfhe.a bin/libthfhe.o

bin/KNN_medical_data: bin/threshold_decryption_functions.o src/KNN_medical_data.cpp
	${TFHE_OPTIONS} g++ -g -fopenmp src/KNN_medical_data.cpp bin/threshold_decryption_functions.o -o bin/KNN_medical_data ${TFHE_LIB} ${BLAS_LIB}

bin/tlwetn: bin/threshold_decryption_functions.o src/TLwe_TN.cpp
	${TFHE_OPTIONS} g++ -g -fopenmp src/TLwe_TN.cpp bin/threshold_decryption_functions.o -o bin/tlwetn ${TFHE_LIB} ${BLAS_LIB}

bin/convert: bin/threshold_decryption_functions.o src/Convert.cpp
	${TFHE_OPTIONS} g++ -g -fopenmp src/Convert.cpp bin/threshold_decryption_functions.o -o bin/convert ${TFHE_LIB} ${BLAS_LIB}

bin/threshold_decryption_functions.o: src/threshold_decryption_functions.cpp
	g++ -g -O3 -c -fopenmp src/threshold_decryption_functions.cpp -o bin/threshold_decryption_functions.o ${BLAS_LIB}

bin/keygen: src/KeyGen.cpp
	${TFHE_OPTIONS} g++ src/KeyGen.cpp -o bin/keygen ${TFHE_LIB}

.PHONY: clean
clean:
	rm -rf bin/* test/*
	touch bin/.gitkeep
	touch test/.gitkeep
	touch test/plain22.txt
	touch test/plain23.txt
	touch test/m0.txt
	touch test/m1.txt
	mkdir -p test/secret.key_shards
	touch test/secret.key_shards/,gitkeep
