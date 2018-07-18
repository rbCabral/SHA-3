all: bench test
CC=clang
keccak_mul_all.o: sources/keccak_mul_all.c
	mkdir -p obj
	$(CC) -c sources/keccak_mul_all.c -march=native  -m64 -O3 -o obj/keccak_mul_all.o

bench.o: sources/bench.c
	mkdir -p obj
	$(CC) -c sources/bench.c -march=native  -m64 -O3 -o obj/bench.o
	
keccak_std.o: sources/sha3_readble.c 
	mkdir -p obj
	$(CC) -c sources/sha3_readble.c -m64 -march=core-avx2 -O -fomit-frame-pointer -o obj/keccak_std.o

bench: 	bench.o keccak_mul_all.o 
	mkdir -p bin
	$(CC) -o bin/SHA-3_bench.x obj/keccak_mul_all.o obj/bench.o -march=native  -m64 -O3 
	
test.o: sources/test.c
	mkdir -p obj
	$(CC) -c sources/test.c -march=native  -m64 -o obj/test.o	
	
test: 	test.o keccak_mul_all.o keccak_std.o
	mkdir -p bin
	$(CC) -o bin/SHA-3_test.x obj/keccak_mul_all.o obj/keccak_std.o obj/test.o -march=native  -m64 	
	
clean:
	rm -rf bin/* obj/*

	
