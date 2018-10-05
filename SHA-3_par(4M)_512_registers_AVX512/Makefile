all: bench test

SHA-3_256.o: sources/SHA-3_256.c
	mkdir -p obj
	gcc -c sources/SHA-3_256.c -march=native  -DKECCAK_256 -m64 -O3 -o obj/SHA-3_256.o

bench.o: sources/bench.c
	mkdir -p obj
	gcc -c sources/bench.c -march=native  -m64 -O3 -o obj/bench.o
	
keccak_std.o: sources/sha3_readble.c 
	mkdir -p obj
	gcc -c sources/sha3_readble.c -m64 -march=core-avx2 -O -fomit-frame-pointer -o obj/keccak_std.o

bench: 	bench.o SHA-3_256.o 
	mkdir -p bin
	gcc -o bin/SHA-3_bench.x obj/SHA-3_256.o obj/bench.o -march=native  -m64 -O3 
	
test.o: sources/test.c
	mkdir -p obj
	gcc -c sources/test.c -march=native  -m64 -o obj/test.o	
	
test: 	test.o SHA-3_256.o keccak_std.o
	mkdir -p bin
	gcc -o bin/SHA-3_test.x obj/SHA-3_256.o obj/keccak_std.o obj/test.o -march=native  -m64 	
	
clean:
	rm -rf bin/* obj/*

	