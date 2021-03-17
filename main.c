/*
*  PROGRAM: pTrimmer
*  FUNCTION: Used to trim off the primer seq from the fastq file
*  INPUT:
*       (1) amplicon file [format: see example]
*       (2) fastq1 and fastq2 [.fq or .fq.gz]
*       (3) output directory [path]
*  OUTPUT:
*       (1) Fastq file after trim off the primer seq [.fq]
*       (2) Amplicon count file [Summary.ampcount]
*
*  AUTHOR: xiaolong Zhang
*  EMAIL:  xiaolongzhang2015@163.com
*  DATE:   2017-09-21 
*  UPDATE: 2021-03-17
*/

#include <pthread.h>
#include <time.h>

#include "utils.h"
#include "query.h"

typedef struct {
    int total, badprim, badqual;
    int ampnum;
    int *ampcount;
} report_t;


uint64_t Total = 0;
uint64_t BadPrim = 0;
uint64_t BadQual = 0;

pthread_mutex_t Lock = PTHREAD_MUTEX_INITIALIZER;


void *process(void *arginfo)
{
    query_t Q;
    report_t rep = {0};
    arginfo_t *Arg = (arginfo_t *)arginfo;
    fastq_t *fq = &Arg->fastq;

    FastqInit(fq, Arg->args, Arg->isread2);
    rep.ampnum = Arg->fwdprim->ampnum;
    err_calloc(rep.ampcount, rep.ampnum, int);

    while (FastqRead(fq)) {
        ++rep.total;
        if (fq->bufnum % BUFNUM == 0) FastqWrite(fq);
        
        Q = PrimQuery(fq->read.seq, Arg);
        if (Q.isfind) /* Even bad quality the read has */
            rep.ampcount[Q.ploc]++;
        else ++rep.badprim;

        Q.badqual = 0; PrimTrim(fq, &Q, Arg);
        if (Q.badqual) rep.badqual++;
    } FastqWrite(fq);

    pthread_mutex_lock(&Lock); 
    amp_t *amp = Arg->fwdprim->amp;
    for (int i=0; i < rep.ampnum; i++)
        amp[i].readnum += rep.ampcount[i];
    
    Total += rep.total; 
    BadPrim += rep.badprim; BadQual += rep.badqual;
    pthread_mutex_unlock(&Lock);

    free(rep.ampcount);
}


int main(int argc, char **argv)
{
    pthread_t tid[2];
    arginfo_t *arginfo[2];
    time_t start, end;
   
    arg_t *args = ParseOpt(argc, argv);
    if (args->help) Usage();

    /* load the primer and kmer index infomation */
    prim_t *fwdprim = GetPrim(args->ampfile);
    prim_t **revprim_list = RevPrim(fwdprim);
    hash_t *fwdindex = InitHash(fwdprim->ampnum << 5);
    PrimIndex(fwdprim, fwdindex, args->kmer);
    hash_t **revindex_list = RevIndex(revprim_list, fwdprim->ampnum, args->kmer);

    /* Sample 1000 reads to calculate the maximum read length and Phred encode */
    int maxrl = ReadLenCheck(args->read1);
    int phred = PhredCheck(args->read1);
    int maxpl = MaxPrimLen(fwdprim);

    for (int i=0; i < args->seqtype+1; ++i) {
        err_calloc(arginfo[i], 1, arginfo_t);
        {
            arginfo[i]->args = args;
            arginfo[i]->fwdprim = fwdprim;
            arginfo[i]->revprim_list = revprim_list;
            arginfo[i]->fwdindex = fwdindex;
            arginfo[i]->revindex_list = revindex_list;
        }
        {
            arginfo[i]->maxpl = maxpl; // maximum primer length
            arginfo[i]->phred = phred; // phred encode mode (phred+33 or phred+64)
            arginfo[i]->maxrl = maxrl; // maximum read length
            arginfo[i]->isread2 = i;
        }
    }

    time(&start);
    
    for (int i=0; i < args->seqtype+1; i++) {
        fprintf(stdout, "[*] Processing the [thread: %d] ...\n", i+1);
        pthread_create(&tid[i], NULL, process, (void*)arginfo[i]);
    }

    for (int i=0; i < args->seqtype+1; i++)
        pthread_join(tid[i], NULL);
    
    /* report the reads count within each amplicon */    
    AmpWrite(fwdprim, args->summary);

    time(&end);
    fprintf(stdout, "Total time consume: %.1f(s)\n", difftime(end,start));

    /* ------------------- summary ------------------- */
    fprintf(stdout, "\n----------------- Summary ------------------------\n");

#ifdef __x86_64__
    fprintf(stdout, "Total reads processed: %lu\n", Total);
    fprintf(stdout, "Reads have bad primer: %lu\n", BadPrim);
    fprintf(stdout, "Reads have bad quality: %lu\n", BadQual);
#else
    fprintf(stdout, "Total reads processed: %llu\n", Total);
    fprintf(stdout, "Reads have bad primer: %llu\n", BadPrim);
    fprintf(stdout, "Reads have bad quality: %llu\n", BadQual);
#endif

    fprintf(stdout, "Reads successfully trimed and have good quality: %.2f %%\n", (float)(Total-BadPrim-BadQual)/Total*100);

    return 0;
}


