# Makefile for pTrimmer
CC = gcc
CFLAGS = -std=c99 -O2 -g
LIBDIR =
LIBS = -lpthread
INCLUDE = 

OBJECT = fastq.o hash.o index.o parse.o query.o dynamic.o main.o

ifeq ($(shell uname -s),Linux)
	PROG = pTrimmer-1.3.1
	LIBS += -lz
else
	PROG = pTrimmer-1.3.1.exe
	INCLUDE += -IWin32
	LIBDIR += -LWin32
	LIBS += -lzdll
endif

$(PROG): $(OBJECT)
	$(CC) $(CFLAGS) $(OBJECT) -o $@ $(INCLUDE) $(LIBDIR) $(LIBS)

fastq.o: fastq.h utils.h
hash.o: hash.h
index.o: hash.h utils.h
parse.o: fastq.h utils.h
query.o: query.h
dynamic.o: dynamic.h
main.o: query.h

.PHONY : clean

clean:
	rm -f $(OBJECT)

