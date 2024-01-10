#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif


#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct {
    size_t l, m;
    char *s;
} kstring_t;
#endif


/* safely open a text file */
#define err_open(_fp, _fn, _mode) do {\
    _fp = fopen(_fn, _mode); \
    if (!_fp) { \
        fprintf(stderr, "\nerr: failed to open %s!\n", _fn); exit(-1); \
    } \
} while(0)

/* safely open a gziped file */
#define err_gzopen(_fp, _fn, _mode) do {\
    _fp = gzopen(_fn, _mode); \
    if (!_fp) { \
        fprintf(stderr, "\nerr: failed to open %s!\n", _fn); exit(-1); \
    } \
} while(0)

/* safely malloc memeory for type */
#define err_malloc(_p, _n, _type) do { \
    _type *tem = (_type *)malloc((_n) * sizeof(_type)); \
    if (!tem) { \
        fprintf(stderr, "\nerr: failed to malloc memory!\n"); exit(-1); \
    } (_p) = tem; \
} while(0)

/* safely calloc memeory for type */
#define err_calloc(_p, _n, _type) do { \
    _type *tem = (_type *)calloc((_n), sizeof(_type)); \
    if (!tem) { \
        fprintf(stderr, "\nerr: failed to calloc memory!\n"); exit(-1); \
    } (_p) = tem; \
} while(0)

/* safely realloc memeory for type */
#define err_realloc(_p, _n, _type) do { \
    _type *tem = (_type *)realloc((_p), (_n) * sizeof(_type)); \
    if (!tem) { \
        fprintf(stderr, "\nerr: failed to realloced memory!\n"); exit(-1); \
    } (_p) = tem; \
} while(0)

#define err_free(_p) do { if (_p) free(_p); } while(0)


/* copy the content of the src to kdest */
kstring_t *k_strcpy(kstring_t *kdest, const char *src);

/* copy n bytes from src to kdest */
kstring_t *k_strncpy(kstring_t *kdest, const char *src, size_t n);


#endif
