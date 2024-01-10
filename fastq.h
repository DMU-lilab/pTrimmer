#ifndef FASTQ_H
#define FASTQ_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include "parse.h"
#include "utils.h"
#include "fileio.h"


#define BUFF_NUM 65536  // maximum number of reads for cache (trimmed reads)


/* fastq type [ read1 or read2 ] */
enum READTYPE { READ1 = 0, READ2 = 1 };

/* illumina quality type phred+33(33<=qual<=75) and phred+64(64<=qual<=106)*/
enum ILLUMINATYPE { Phrd33 = 0, Phrd64 = 1};


/*! @typedef read_t
  @abstract the new read object
  @field  name      the name of read('@' will be repleaced with pair_marker if 'pair_check' is given).
                    eg. @ST-E00126:HWM7:3173:1784 2:N:AAGAC -> 2ST-E00126:HWM7:3173:1784
  @field  seq       the sequence of the read
  @field  comment   the comment of the read, usually is a character of '+'
  @field  qual      the quality of the read
 */
typedef struct __read_t {
    kstring_t name;
    kstring_t seq;
    kstring_t comment;
    kstring_t qual;
} read_t;


/*! @typedef fastq_t
 @abstract structure for the single fastq file
 @field in_hd        the input stream handle of the input fastq file
 @field out_hd       the output stream handle of the trimmed fastq file
 @field read         a fastq read object
 @field size         the current number of reads in the cache
 @field max          the memory allocated for the cache
 @field cache        the buffer for the fastq cache
*/
typedef struct __fastq_t {
    GzStream *in_hd;
    GzStream *out_hd;
    read_t read;
    int size, max;
    read_t *cache;
} fastq_t;


/* prototype function */
void FastqInit(fastq_t *, arg_t *, int);
void FastqDestroy(fastq_t *);
int PhredCheck(char *);
int ReadLenCheck(char *);
float MeanQuality(kstring_t *, int);
int FastqGetRead(GzStream *, read_t *);
void FastqWriteCache(fastq_t *);
void FastqFreeRead(read_t *);

/* parse the comand line parmerters*/
void Usage(void);
arg_t *ParseOpt(int, char **);

#endif
