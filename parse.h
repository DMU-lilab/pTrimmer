/*************************************************************************
    > File Name: parse.h
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com
    > Created Time: 2024年01月10日 星期三 10时57分32秒
 ************************************************************************/

#ifndef PTRIMMER_PARSE_H
#define PTRIMMER_PARSE_H

#include <getopt.h>
#include <string.h>

#ifdef _WIN32
  #define PATH_MAX 260
#elif __linux__
  #include <linux/limits.h>
#else
#include <zlib.h>
  #include <limits.h>
#endif

/* sequencing type SE (Single-end) and PE (Paired-end) */
enum SEQTYPE { SE = 0, PE = 1 };


/*! @typedef arg_t
 @abstract structure for the comand line args
 @field help         [0 | 1] if 1: print the help infomation
 @field kmer         the length for primer index
 @field mismatch     mismatch allowed between the primer and sequence
 @field seqtype      sequencing type which include single-end and paired-end
 @field keep         if true, keep the complete reads that can not find primer
 @field gzip         if given, compress the trimmed reads in Gzip format
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
    int gzip;  // 0 (not gzipped) or 1 (gzipped)
    int minqual; // default: 20
    char ampfile[PATH_MAX];
    char read1[PATH_MAX];
    char read2[PATH_MAX];
    char trim1[PATH_MAX];
    char trim2[PATH_MAX];
    char summary[PATH_MAX];
} arg_t;


#endif //PTRIMMER_PARSE_H
