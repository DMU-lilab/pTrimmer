/* seq query function */

#include "query.h"
#include "utils.h"
#include "dynamic.h"

/* buffer length for kmer search */
#define BLEN 6

/*! @funciton: check the mismatch base
 *   @parmeters:
 *   @    str1      the pointer to the string 1 [char *]
 *   @    str2      the pointer to the string 2 [char *]
 *   @    mismatch  the maxmum mismatch base allowed [int]
 *   @return:
 *   @    0     False
 *   @    1     True
*/
static int MisCheck(char *str1, char *str2, int mismatch)
{
    int misnum = 0;
    
    if (!*str1 || !*str2) return 0;

    while (*str1 && *str2) {
        if (*str1++ != *str2++)
            misnum++;

        if (misnum > mismatch)
            return 0;
    }
    return 1;
}

/*! @funciton: calculate the maximux primer length
 *   @parmeters:
 *   @    P        pointer to the primer struct [prim_t *]
 *   @return:
 *   @    max      the maximum primer length
*/
int MaxPrimLen(prim_t *P)
{
    int max=0, flen, rlen;
    amp_t *a;

    for (int i=0; i < P->ampnum; ++i) {
        a = &P->amp[i];
        flen = (int)strlen(a->fwdprim);
        rlen = (int)strlen(a->revprim);
        max = flen > max ? flen : max;
        max = rlen > max ? rlen : max;
    }

    if (max > MAX_PRIMMER_LEN) {
        fprintf(stderr, "[Err:%s] the primer length can not longer than %d!", __func__, max);
        exit(-1);
    }

    return max;
}


query_t HammingMode(char *seq, prim_t *P, int mis, int index, status *S)
{
    query_t Q = {0};

    for (int i=0; i < S->node->num; ++i) {
        loc_t *loc = &S->node->loc[i];
        int shift = index - loc->sloc;

        char *string = seq;
        char *pstr = loc->frloc ?
            P->amp[loc->ploc].revprim : P->amp[loc->ploc].fwdprim;

        if (shift >= 0) {
            /* string: ATGCCGAATTGATGCCTGA
             *   pstr:   GCCGAA
             * index=2; sloc=0 => shift=2 */
            string = seq + shift;
        }
        else if (shift >= -mis) {
            /* string:   GCCGAATTGATGCCTGA
             *   pstr: ATGCCGAA
             * index=0; sloc=2 => shift=-2 */
            pstr -= shift;
        }
        else continue;

        if (MisCheck(string, pstr, mis)) {
            Q.isfind = 1;
            Q.ploc = loc->ploc;
            Q.pstart = shift < 0 ? 0 : shift;
            Q.pend = Q.pstart + (int)strlen(pstr) - 1;
            return Q;
        }
    }
    return Q;
}

query_t DynamicMode(char *seq, compo_t *compo, prim_t *P, int mis)
{
    dynamic_t d;
    query_t Q = {0};

    char *pstr = compo->frloc ?
            P->amp[compo->ploc].revprim : P->amp[compo->ploc].fwdprim;

    int slen = (int)strlen(seq);
    int plen = (int)strlen(pstr);
    /* search with dynamic algorithm */
    d = global_align(seq, slen, pstr, plen, mis);

    if (d.matches >= plen-mis && d.errors <= mis) {
        Q.isfind = 1;
        Q.ploc = compo->ploc;
        Q.pstart = d.start1; Q.pend = d.stop1 -1;
    }
    return Q;
}

static int _compscore(const void *a, const void *b) {
    /* sorted by descending order */
    return ((compo_t*)b)->score - ((compo_t*)a)->score;
}

static void CalHits(hit_t *hit, node_t *node, int seqi)
{
    loc_t *loc;
    compo_t *c;

    for (int i=0; i < node->num; ++i) {
        loc = &node->loc[i];

        if (hit->hitnum % 0x8 == 0) {
            /* when n==0, realloc is similar to malloc */
            int n = hit->hitnum + 0x8;
            err_realloc(hit->compo, n, compo_t);
        }

        c = hit->compo;
        int flag = 0;
        for (int j=0; j < hit->hitnum; ++j) {
            if (loc->ploc == c[j].ploc && loc->frloc == c[j].frloc) {
                /* the loc has existed in the hit */
                c[j].score++; 
                flag = 1; break;
            }
        }

        if (!flag) { /* the loc not existed in the hit */
            c = &hit->compo[hit->hitnum++];
            c->ploc = loc->ploc; c->frloc = loc->frloc;
            c->sloc = loc->sloc; c->seqi = seqi;
            c->score = 1; /* the value may not be zero after relloc */
        }
    }
}


static hit_t *GetHits(char *seq, hash_t *H, int kmer)
{
    char key[KEYLEN];
    status S;
    hit_t *hit;

    err_calloc(hit, 1, hit_t);
    int cutlen = (int)strlen(seq);

    /* 'N' is imposible occured in primer index */
    for (int i=0; i < cutlen; ++i) {
        if (seq[i] == 'N') continue;

        strncpy(key, seq+i, kmer); key[kmer] = '\0';
        Search(key, H, &S);
        if (S.find == 0) continue;

        CalHits(hit, S.node, i);
    }

    if (hit->hitnum)
        /* sort the score by descending order */
        qsort(hit->compo, hit->hitnum, sizeof(compo_t), _compscore);

    return hit;
}


/* CORE FUNCTION 
 *      seq: AGAAATTTGCGGAGTAAGTTGCGCTGGGGCTTTCGGCGGCGGCGATTT
 *   primer:      TTTGCGGA
 *                |      |
 *               pstart  pend
 * */
static query_t SeqQuery(char *seq, hash_t *H, prim_t *P, int mis, int kmer)
{
    status S;
    query_t Q = {0};
    char key[KEYLEN];

    for (int i=0; i < (BLEN<<1); i++) {
        /* 'N' is imposible occured in primer index */
        if (*(seq+i) == 'N') continue;

        strncpy(key, seq+i, kmer); key[kmer] = '\0';
        Search(key, H, &S);
        if (S.find == 0) continue;

        /* find the primer seq with Hanming Mode */
        Q = HammingMode(seq, P, mis, i, &S);
        if (Q.isfind) return Q;
    }

    /* Start the Dynamic Mode */
    hit_t *hit = GetHits(seq, H, kmer);
    if (!hit->hitnum) { /* can't find even one kmer hit */
        err_free(hit); return Q;
    }
    Q = DynamicMode(seq, &hit->compo[0], P, mis);
    err_free(hit->compo); err_free(hit);

    return Q;
}


/*! @funciton: get the insert fragment start and end (0-base)
 *   @parmeters:
 *   @    seq       the pointer to fastq sequence [char *]
 *   @    arg       all the necessary parmeters [arginfo_t *]
 *   @return:
 *   @    Q         the start and end of the insert fragment [query_t]
*/
query_t PrimQuery(char *seq, arginfo_t *arg)
{
    query_t f, r, Q ;
    char temseq[MAX_PRIMMER_LEN<<1];

    strncpy(temseq, seq, arg->maxpl+BLEN);
    temseq[arg->maxpl+BLEN] = '\0';
    f = SeqQuery(temseq, arg->fwdindex, 
                arg->fwdprim, arg->args->mismatch, arg->args->kmer);
    if (!f.isfind) {
        Q.isfind = 0; return Q; 
    }

    amp_t *amp = &(arg->fwdprim->amp[f.ploc]);
    if (f.pend+amp->insertlen > arg->maxrl) {
        Q.isfind = 1;
        Q.ploc = f.ploc;
        Q.pstart = f.pend + 1;
        Q.pend = 0; // the revprimer is not sequenced!
        return Q;
    }

    int rstart = f.pend + amp->insertlen - BLEN;
    r = SeqQuery(seq+rstart, arg->revindex_list[f.ploc],
                arg->revprim_list[f.ploc], arg->args->mismatch, arg->args->kmer);
    if (!r.isfind) { 
        Q.isfind = 0; return Q; 
    }

    Q.isfind = 1;
    Q.ploc = f.ploc;
    Q.pstart = f.pend + 1;
    Q.pend = rstart + r.pstart -1;

    return Q;
}

/* add primer information to comment of the trimmed read */
void AddPrimerInfo(kstring_t *comment, query_t *Q, amp_t *amp)
{
    char buf[64];
    int f_start=0, f_end=0, r_start=0, r_end=0;

    /* failed to find the read's primer */
    if (Q->isfind == 0) {
        sprintf(buf, "+ P=%d F=%d:%d R=%d:%d\n", 0, f_start, f_end, r_start, r_end);
        k_strcpy(comment, buf);
        return ;
    }

    /* calculate the primer location on the read */
    amp_t *amp_obj = &amp[Q->ploc];

    f_start = Q->pstart - amp_obj->fwd_len + 1;
    f_end = Q->pstart;

    if (Q->pend > 0) {  /* read-through condition */
        r_start = Q->pend + 2;
        r_end = r_start + amp_obj->rev_len - 1;
    }
    sprintf(buf, "+ P=%d F=%d:%d R=%d:%d\n", Q->ploc+1, f_start, f_end, r_start, r_end);
    k_strcpy(comment, buf);
}


/*! @funciton: trim off the primer seq
 *   @parmeters:
 *   @    fq        the pointer to fastq structure [fastq_t *]
 *   @    Q         the query result from PrimQuery [query_t *]
 *   @    arg       the args required [arginfo_t *]
 *   @return:
 *   @    void
*/
void PrimTrim(fastq_t *fq, query_t *Q, arginfo_t *arg)
{
    read_t *dest_read, *src_read;
    amp_t *amp_obj = arg->fwdprim->amp;
    float min_qual = (float)arg->args->minqual;

    dest_read = &fq->cache[fq->size++];
    src_read = &fq->read;
    k_strcpy(&dest_read->name, src_read->name.s);

    /* add primer info for trimmed read in field of comment */
    if (arg->args->info) AddPrimerInfo(&dest_read->comment, Q, amp_obj);
    else k_strcpy(&dest_read->comment, "+\n");

    /* trimming the primer sequence */
    while (Q->isfind) {  // find the primer sequence
        k_strcpy(&dest_read->seq, src_read->seq.s+Q->pstart);
        k_strcpy(&dest_read->qual, src_read->qual.s+Q->pstart);
        if (Q->pend) {
            /* read through */
            int inlen = Q->pend - Q->pstart +1;
            k_strncpy(&dest_read->seq, dest_read->seq.s, inlen);
            k_strcat(&dest_read->seq, "\n");
            k_strncpy(&dest_read->qual, dest_read->qual.s, inlen);
            k_strcat(&dest_read->qual, "\n");
        }
        if (dest_read->seq.l < 1) break;  /* the length is too short */

        /* evaluate the quality */
        if (MeanQuality(&dest_read->qual, arg->phred) < min_qual) {
            Q->badqual = 1; break;
        }
        return ;
    }

    /* failed to find the primer sequence or bad quality or too short length */
    if (arg->args->keep) {
        k_strcpy(&dest_read->seq, src_read->seq.s);
        k_strcpy(&dest_read->qual, src_read->qual.s);
        return ;
    }

    /* replace the sequence and quality */
    k_strcpy(&dest_read->seq, "NNNNNNNNNNNNNNNNNNNN\n");
    k_strcpy(&dest_read->qual, "!!!!!!!!!!!!!!!!!!!!\n");
}

