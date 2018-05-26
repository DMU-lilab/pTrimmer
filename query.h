#ifndef QUERY_H
#define QUERY_H

#include "index.h"
#include "fastq.h"

/*! @typedef compo_t
 @abstract structure for the kmer located at same primer sequence
 @field ploc        the primer index in the amplicon list
 @field frloc       forward or revrse index [0=>forward and 1=>reverse]
 @field sloc        the index of key in the primer sequence
 @field seqi        the index of 1-th kmer occoured in the seq
 @field score       number of kmer from this primer sequence
*/
typedef struct {
    int ploc, frloc, sloc;
    int seqi;
    int score;
} compo_t;

/*! @typedef hit_t
 @abstract structure for all the possible primer hit
 @field hitnum      number of potential primer hit
 @field compo       the array for each primer kmer hit
*/
typedef struct __hit_t {
    int hitnum;
    compo_t *compo;
} hit_t;

/*! @typedef query_t
 @abstract structure for the return value of query function
 @field isfind      whether find the match primer
 @field badqual     whether the reads is bad quality after trimming
 @field ploc        the primer index in the amplicon list
 @field pstart      primer start index in the query sequecing
 @field pend        primer end index in the query sequecing
*/
typedef struct __query_t {
    int isfind;
    int badqual;
    int ploc;
    int pstart;
    int pend;
} query_t;


/*! @typedef arginfo_t
 @abstract structure for the thread parmeters
 @field isread2         1=> read2 otherwise 0=> read1
 @field maxrl           maximum read length (sample from 1000 reads)
 @field maxpl           maximum primer length
 @field args            come from comand line parmeters
 @field fwdindex        foroward primer index
 @field revindex_list   reverse primer index table list
 @field fwdprim         foroward primer list
 @field revprim_list    reverse primer list
 @field fastq           read1 or read2 fastq container
*/
typedef struct __arginfo_t {
    int isread2;
    int maxrl, maxpl;
    int phred;
    arg_t *args;
    hash_t *fwdindex;
    hash_t **revindex_list;
    prim_t *fwdprim;
    prim_t **revprim_list;
    fastq_t fastq;
} arginfo_t;


/* prototype function */
int MaxPrimLen(prim_t *);
query_t PrimQuery(char *, arginfo_t *);
void PrimTrim(fastq_t *, query_t *, arginfo_t *);

#endif
