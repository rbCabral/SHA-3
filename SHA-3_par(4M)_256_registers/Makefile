all: bench test

keccak_mul_all.o: sources/keccak_mul_all.c
	gcc -c sources/keccak_mul_all.c -march=native  -m64 -O3 -o obj/keccak_mul_all.o

bench.o: sources/bench.c
	gcc -c sources/bench.c -march=native  -m64 -O3 -o obj/bench.o
	
keccak_std.o: sources/sha3_readble.c 
	gcc -c sources/sha3_readble.c -m64 -march=core-avx2 -O -fomit-frame-pointer -o obj/keccak_std.o

bench: 	bench.o keccak_mul_all.o 
	gcc -o bin/SHA-3_bench.x obj/keccak_mul_all.o obj/bench.o -march=native  -m64 -O3 
	
test.o: sources/test.c
	gcc -c sources/test.c -march=native  -m64 -o obj/test.o	
	
test: 	test.o keccak_mul_all.o keccak_std.o
	gcc -o bin/SHA-3_test.x obj/keccak_mul_all.o obj/keccak_std.o obj/test.o -march=native  -m64 	
	
clean:
	rm -rf bin/* obj/*

	