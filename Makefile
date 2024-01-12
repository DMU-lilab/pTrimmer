# Makefile for pTrimmer
# update: 20240110

DEBUG = 0

CC = gcc
CFLAGS = -std=c99
LIBDIR =
LIBS = -lpthread -lz
INCLUDE = 

ifeq ($(DEBUG), 1)
    CFLAGS += -g -O0 # enable debugging
else
    CFLAGS += -O2
endif


OBJECT = fastq.o hash.o index.o parse.o query.o dynamic.o fileio.o utils.o main.o

ifeq ($(OS),Windows_NT)
	PROG = pTrimmer.exe
	LIBS += -static
else  # Linux and Darwin
	PROG = pTrimmer
endif

all: $(PROG)

$(PROG): $(OBJECT)
	$(CC) $(CFLAGS) $(OBJECT) -o $@ $(INCLUDE) $(LIBDIR) $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

.PHONY : clean

clean:
	$(RM) -f $(OBJECT) $(PROG)

