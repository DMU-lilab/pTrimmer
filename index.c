/* building index with hash table */

#include "index.h"
#include "utils.h"

#ifdef _WIN32
  #define PATH_MAX 260
  #define _PDELIM_ 92 // '\'
#elif __linux__
  #include <linux/limits.h>
  #define _PDELIM_ 47 // '/'
#else
  #include <limits.h>
  #define _PDELIM_ 47 // '/'
#endif


/*! @funciton: split the primer line
 *   @parmeters:
 *   @    filename   the pointer to the filename of the primer [char *]
 *   @return:
 *   @    p          the pointer to the primer structure [prim_t *]
*/
prim_t *GetPrim(char *filename)
{
    char buf[BUFLINE];
    amp_t *a; prim_t *p;
    FILE *file;

    err_open(file, filename, "r");
    err_calloc(p, 1, prim_t);
    err_calloc(p->amp, 1, amp_t);

    while (fgets(buf, BUFLINE, file)) {
        if (buf[0] != '#') {
            /* realloc necessary memory for amplicon list */
            if (p->ampnum % PRIMNUM == 0) {
                int newsize = p->ampnum + PRIMNUM;
                err_realloc(p->amp, newsize, amp_t);
            }
            a = &(p->amp[p->ampnum]);
            a->readnum = 0;
            sscanf(buf, "%s%s%d%s", a->fwdprim, a->revprim, &a->insertlen, a->auxinfo);
            p->ampnum++;
        }
    } fclose(file);

    return p;
}


/*! @funciton: get the suffix of path
 *   @parmeters:
 *   @    path       the pointer to the path [char *]
 *   @return:
 *   @               the pointer to the suffix start
 *   NOTE: Deprecated Now!!!
*/
static char *GetSuffix(char *path)
{
    int plen, i;

    plen = strlen(path);
    for (i=plen; i > 0; i--) {
        /* get base path for out amplicon count */
        /*　_t => xlzh_trim_fq */
        if (path[i] == '_' && path[i+1] == 't')
            break;
    }
    return (path+i);
}


/*! @funciton: get the reverse and complement sequecne
 *   @parmeters:
 *   @    string     the pointer to target string [char *]
 *   @return:
 *   @    void       the pointer to the suffix start
*/
void RevComp(char *string)
{
    uint8_t tmp;
    char *last;

    for (last=string; *last; ++last) ;
    string--;
    while (++string <= --last) {
        tmp = *string;
        *string = _BASE[*last];
        *last = _BASE[tmp];
    }
    return ;
}


/*! @funciton: get the reverse primer structue
 *   @parmeters:
 *   @    P             the pointer to forward primer [prime_t *]
 *   @return:
 *   @   primlist       the pointer to the reverse primer list [prim_t **]
*/
prim_t **RevPrim(prim_t *P)
{
    prim_t **primlist;
    
    err_calloc(primlist, P->ampnum, prim_t*);
    for (int i=0; i < P->ampnum; i++) {
        err_malloc(primlist[i], 1, prim_t);
        err_malloc(primlist[i]->amp, 1, amp_t);

        primlist[i]->ampnum = 1;
        primlist[i]->amp[0] = P->amp[i];
        RevComp(primlist[i]->amp[0].fwdprim);
        RevComp(primlist[i]->amp[0].revprim);
    }

    return primlist;
}


/*! @funciton: build index for the primer structure
 *   @parmeters:
 *   @    p          the pointer to primer [prime_t *]
 *   @    T          the hast talbe contain the index [hash_t *]
 *   @    kmer       the index length [hash_t *]
 *   @return:
 *   @    void
*/
void PrimIndex(prim_t *p, hash_t *T, int kmer)
{
    loc_t loc;
    char slic[KEYLEN];

    for (int i=0; i < p->ampnum; i++) {
        amp_t *amp = &p->amp[i];

        for (int j=0; j < 2; j++) {
            char *primseq = j ? amp->revprim : amp->fwdprim;
            int plen = strlen(primseq) -kmer + 1;

            for (int k=0; k < plen; k++) {
                strncpy(slic, primseq+k, kmer);
                slic[kmer] = '\0';
                loc.ploc = i; loc.frloc = j; loc.sloc = k;
                Insert(slic, &loc, T);
            }
        }   
    }
    return ;
}


/*! @funciton: get the index list of amplicon reverse primer
 *   @parmeters:
 *   @    P             the pointer to the primer list [prim_t **]
 *   @    ampnum        the amplicon count [int]
 *   @    kmer          the kmer for indexing [int]
 *   @return:
 *   @    revlist       [hash_t **]
*/
hash_t **RevIndex(prim_t **P, int ampnum, int kmer)
{
    hash_t **revlist;

    err_calloc(revlist, ampnum, hash_t*);
    for (int i=0; i < ampnum; i++) {
        revlist[i] = InitHash(1<<5);
        PrimIndex(P[i], revlist[i], kmer);
    }

    return revlist;
}


/*! @funciton: write out the amplicon infomation
 *   @parmeters:
 *   @    prim          the pointer to the forward primer [prim_t *]
 *   @    path          the pointer to the summary file path [char *]
 *   @return:
 *   @    void        
*/
void AmpWrite(prim_t *prim, char *path)
{
    FILE *fp;

    if (!path[0]) // summary path is not passed by user
        err_open(fp, "Summary.ampcount", "w");
    else
        err_open(fp, path, "w");
   
    fprintf(fp, "FwdPrimer\tRevPrimer\tAmpCount\tAuxInfo\n");
    for (int i=0; i < prim->ampnum; i++) {
        amp_t *amp = &prim->amp[i];
        fprintf(fp, "%s\t%s\t%d\t%s\n", \
                amp->fwdprim, amp->revprim, amp->readnum, amp->auxinfo);
    } fclose(fp);

    return ;
}
