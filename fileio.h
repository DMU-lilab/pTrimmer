/*************************************************************************
    > File Name: fileio.h
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com
    > Created Time: 2024年01月10日 星期三 12时10分50秒
 ************************************************************************/

#ifndef PTRIMMER_FILEIO_H
#define PTRIMMER_FILEIO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"  // for kstring_t

#ifdef _WIN32
#include "Win32\zlib.h"
#elif __linux__
#include <zlib.h>
#else
#include <zlib.h>
#endif


/* number of bytes cached from *.gz file */
#define GZ_BUFF_SIZE 1048576
#define READ_LEN_MAX 1048576


/*! @typedef GzStream
 @abstract structure for gz/bz2 handle
 @field bz2_fp       file handle of *.bgz2
 @field gz_fp        file handle of *.gz
 @field out_fp       file handle of output file
 @field buf          used to store decompressed data, whose size is 'gz_stream_buff_size'
 @field begin, end   begin and end index in the buf
 @field is_write     is_write:1 -> (for file write)'w', otherwise -> 'r' (for file read)
 @field is_eof       is_eof:1 -> the end of the file
*/
typedef struct {
    gzFile gz_fp;
    FILE *out_fp;
    char *buf;
    int bzerror;
    int begin;
    int end;
    int is_write;
    int is_eof;
} GzStream;


/*! @function: open a gz or bz2 file
  @param  file        input file name
  @param  mode        operation can only be 'r'(read) and 'w'(write)
  @return             GzStream object
 */
GzStream *gz_stream_open(char *file, char *mode);


/*! @function: read one line(until delimiter) from the compressed file
  @param  gz          GzStream object
  @param  delimiter   delimiter used to decide where to break, e.g. '\n'
  @param  ks_str      kstring_t type of string (must be NULL for the first time)
  @return             operation status: 0->EOF; 1->OK; -1:truncate
 */
int gz_read_util(GzStream *gz, char delimiter, kstring_t *ks_str);


/*! @function:  write one line to the gzip compress file
  @param  gz          GzStream object
  @param  ks_str      kstring_t type of string (must be NULL for the first time)
  @return             void
 */
void gz_write_util(GzStream *gz, kstring_t *ks_str);


/*! @function: destroy the GzStream object
  @param  gz          GzStream object
  @return             void
 */
void gz_stream_destroy(GzStream *gz);



#endif //PTRIMMER_FILEIO_H
