# Makefile for pTrimmer
# update: 20240110

DEBUG = 0

CC = gcc
CFLAGS = -std=c99
LIBDIR =
LIBS = -lpthread
INCLUDE = 

ifeq ($(DEBUG), 1)
    CFLAGS += -g -O0 # enable debugging
else
    CFLAGS += -O2
endif


OBJECT = fastq.o hash.o index.o parse.o query.o dynamic.o main.o

ifeq ($(shell uname -s),Linux)
	PROG = pTrimmer
	LIBS += -lz
	RM = rm
else ifeq ($(shell uname -s),Darwin)
	PROG = pTrimmer
	LIBS += -lz
	RM = rm
else
	PROG = pTrimmer.exe
	INCLUDE += -IWin32
	LIBDIR += -LWin32
	LIBS += -lzdll
	RM = del
endif

all: $(PROG)

$(PROG): $(OBJECT)
	$(CC) $(CFLAGS) $(OBJECT) -o $@ $(INCLUDE) $(LIBDIR) $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

.PHONY : clean

clean:
	$(RM) -f $(OBJECT) $(PROG)

