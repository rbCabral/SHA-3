all: bench test
CC=clang

SHA-3_256.o: sources/SHA-3_256.c
	mkdir -p obj
	$(CC) -c sources/SHA-3_256.c -march=native  -m64 -O3 -o obj/SHA-3_256.o -g -Wall -Wextra

bench.o: sources/bench.c
	mkdir -p obj
	$(CC) -c sources/bench.c -march=native  -m64 -O3 -o obj/bench.o -g -Wall -Wextra
	
keccak_std.o: sources/sha3_readble.c 
	mkdir -p obj
	$(CC) -c sources/sha3_readble.c -m64 -march=core-avx2 -O -fomit-frame-pointer -o obj/keccak_std.o -g -Wall -Wextra

bench: 	bench.o SHA-3_256.o 
	mkdir -p bin
	$(CC) -o bin/SHA-3_bench.x obj/SHA-3_256.o obj/bench.o -march=native  -m64 -O3 -g -Wall -Wextra
	
test.o: sources/test.c
	mkdir -p obj
	$(CC) -c sources/test.c -march=native  -m64 -o obj/test.o -g -Wall -Wextra
	
test: 	test.o SHA-3_256.o keccak_std.o
	mkdir -p bin
	$(CC) -o bin/SHA-3_test.x obj/SHA-3_256.o obj/keccak_std.o obj/test.o -march=native  -m64 -g -Wall -Wextra

runT: test
	bin/SHA-3_test.x 4096

runB: bench
	bin/SHA-3_bench.x 4096

clean:
	rm -rf bin/* obj/*

	
