CC?=gcc

all: crimp

crimp: src/crimp.c
	$(CC) -std=c99 -Wall -O3 -o crimp src/crimp.c -lm
