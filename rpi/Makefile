CC = gcc

TP ?= 47

DEPS = build/lwe.o build/rlwe.o build/random.o build/perf.o

all: main

main: src/main.c $(DEPS)
	$(CC) -g -Iinclude -DWITHOS -DTP=$(TP) src/main.c $(DEPS) -o main -lc -lm

build/lwe.o: src/tfhe/lwe.c
	mkdir -p build
	$(CC) -c -g -Iinclude src/tfhe/lwe.c -o build/lwe.o -lc -lm

build/rlwe.o: src/tfhe/rlwe.c
	mkdir -p build
	$(CC) -c -g -Iinclude src/tfhe/rlwe.c -o build/rlwe.o -lc -lm

build/random.o: src/utils/random.c
	mkdir -p build
	$(CC) -c -g -Iinclude src/utils/random.c -o build/random.o -lc -lm

build/perf.o: src/utils/perf.S
	mkdir -p build
	$(CC) -c -g -Iinclude src/utils/perf.S -o build/perf.o -lc -lm

.PHONY: clean
clean:	
	rm -rf main build