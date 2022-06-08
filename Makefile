CC=gcc
RM=rm
MAKE=make
CFLAGS=-std=c99 -O3 -pthread
all:
	$(MAKE) --no-print-directory ./bin/crisflash

./bin/crisflash: ./src/main_crisflash.c ./src/read.c ./src/nary_tree.c ./src/vcf.c
	mkdir -p bin
	$(CC) $(CFLAGS) ./src/main_crisflash.c ./src/read.c ./src/nary_tree.c ./src/vcf.c -o ./bin/crisflash

clean: ./bin/crisflash
	$(RM) ./bin/crisflash
