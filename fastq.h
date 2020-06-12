#ifndef FASTQ_H
#define FASTQ_H

#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>

#ifdef _WIN32
  #include "Win32\zlib.h"
  #define PATH_MAX 260
  #define BUFNUM 1024 // default buf number for the fastq cache
#elif __linux__
  #include <zlib.h>
  #include <linux/limits.h>
  #define BUFNUM 2048
#else
  #include <zlib.h>
  #include <limits.h>
  #define BUFNUM 2048
#endif


#define FQLINE 256    // maxmum length for the fastq line
#define BUFLINE 1024  // maxmum bufer line for file reading


/* fastq type [ read1 or read2 ] */
enum READTYPE { READ1 = 0, READ2 = 1 };
/* sequencing type SE (Single-end) and PE (Paired-end) */
enum SEQTYPE { SE = 0, PE = 1 };
/* illumina quality type phred+33(33<=qual<=75) and phred+64(64<=qual<=106)*/
enum ILLUMINATYPE { Phrd33 = 0, Phrd64 = 1};

/*! @typedef arg_t
 @abstract structure for the comand line args
 @field help         [0 | 1] if 1: print the help infomation
 @field kmer         the length for primer index
 @field mismatch     mismatch allowed between the primer and sequence
 @field seqtype      sequencing type which include single-end and paired-end
 @field keep         if true, keep the complete reads that can not find primer
 @field minqual      the mimimum average quality of the reads after trimming
 @field ampfile      the path of amplicon file [format: see example]
 @field read1        the path of fastq file of R1
 @field read2        the path of fastq file of R2
 @field trime1       the path of the trimed fastq file of R1
 @field trime2       the path of the trimed fastq file of R2
 @field summary      the path of the summary file [default: Summary.ampcount]
*/
typedef struct __arg_t {
    int help;
    int kmer;
    int mismatch;
    int seqtype; // SE(single) or PE(pair)
    int keep;    // 0 or 1
    int minqual; // default: 20
    char ampfile[PATH_MAX];
    char read1[PATH_MAX];
    char read2[PATH_MAX];
    char trim1[PATH_MAX];
    char trim2[PATH_MAX];
    char summary[PATH_MAX];
} arg_t;


/*! @typedef read_t
 @abstract structure for one fastq read group
 @field name         the seqname for the read
 @field seq          the sequence for the read
 @field mark         the fixed marker['+']
 @field qual         the quality for the read
*/
typedef struct __read_t {
    char name[FQLINE];
    char seq[FQLINE];
    char mark[FQLINE];
    char qual[FQLINE];
} read_t;


/*! @typedef fastq_t
 @abstract structure for the single fastq file
 @field in, out      the input and output pointer for the fastq [FILE *]
 @field read         a fastq read group
 @field bufnum       the bufnum for the cache
 @field cache        the buffer for the fastq cache
*/
typedef struct __fastq_t {
    gzFile in;
    FILE *out;
    read_t read;
    int bufnum;
    read_t cache[BUFNUM];
} fastq_t;


/* prototype function */
void FastqInit(fastq_t *, arg_t *, int);
int PhredCheck(char *);
int ReadLenCheck(char *);
float MeanQuality(char *, int);
int FastqRead(fastq_t *);
void FastqWrite(fastq_t *);

/* parse the comand line parmerters*/
void Usage(void);
arg_t *ParseOpt(int, char **);

#endif
