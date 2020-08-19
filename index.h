#ifndef INDEX_H
#define INDEX_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdlib.h>

#include "hash.h"

#define MAXPRIMLEN 128 // ideal maxmum length for the primer seq
#define PRIMNUM 128    // initial primer number
#define AUXLEN 512     // maxmum length for auxinfo
#ifndef BUFLINE
    #define BUFLINE 1024
#endif

/* A:T; T:A, G:C, C:G, N:N */
static  uint8_t _BASE[96] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,84,0,71,0,0,0,67,0,0,0,0,0,0,78,0,
    0,0,0,0,65,0,0,0,0,0,0,0,0,0,0,0
};


/*! @typedef amp_t
 @abstract structure for the amplicon line
 @field insertlen      the insert length for the paired primer
 @field readnum        the read number within the amplicon
 @field fwdprim        forward primer sequence
 @field revprim        reverse primer sequence
 @field auxinfo        other filed infomation
*/
typedef struct __amp_t {
    int insertlen;
    int readnum;
    char fwdprim[MAXPRIMLEN];
    char revprim[MAXPRIMLEN];
    char auxinfo[AUXLEN];
} amp_t;


/*! @typedef prim_t
 @abstract structure for the primer
 @field ampnum     the amplicon number with the primer structure
 @field maxplen    the maximum primer length
 @field amp        the pointer to the amplicon list
*/
typedef struct __prim_t {
    int ampnum;
    int maxplen;
    amp_t *amp;
} prim_t;


/* prototype function */
prim_t *GetPrim(char *);
prim_t **RevPrim(prim_t *);
hash_t **RevIndex(prim_t **, int, int);
void PrimIndex(prim_t *, hash_t *, int);
void AmpWrite(prim_t *, char *);

#endif

