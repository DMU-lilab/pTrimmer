#include "fastq.h"
#include "utils.h"

#ifdef _WIN32
  #define _PDELIM_ 92 // '\'
#else
  #define _PDELIM_ 47 // '/'
#endif

/*! @funciton: get the basename of path
 *   @parmeters:
 *   @    path       the pointer to the path [char *]
 *   @return:
 *   @               the pointer to the basename start
*/
static char *BaseName(char *path)
{
    int plen, i;

    plen = strlen(path);
    for (i=plen; i > 0; i--) {
        if (path[i] == _PDELIM_)
            break;
    }
    if (i) return (path+i+1);
    else return path;
}


/*! @funciton: fastq initialization
 *   @parmeters:
 *   @    fq         the pointer to the fastq_t structure [fastq_t *]
 *   @    args       the pointer to the arg_t structure [fastq_t *]
 *   @    type       the read1 or read2 for the fastq file
 *   @return:
 *   @    void
*/
void FastqInit(fastq_t *fq, arg_t *args, int type)
{
    char bname[128];
    char outname[PATH_MAX];

    fq->bufnum = 0;

    if (type == READ1) {
        err_gzopen(fq->in, args->read1, "r");
        err_open(fq->out, args->trim1, "w");
    }
    else {
        err_gzopen(fq->in, args->read2, "r");
        err_open(fq->out, args->trim2, "w");
    }
    gzbuffer(fq->in, 524288);
}


/*! @funciton: Phred encoding check
 *   @parmeters:
 *   @    infile     the pointer to the input fastq file [char *]
 *   @return:
 *   @    Phred33[0] or Phred64[1]
*/
int PhredCheck(char *infile)
{
    char *b;
    fastq_t fq;
    int min=128, max=0;

    err_gzopen(fq.in, infile, "r");
    for (int i=0; i < 1000; ++i) {
        if (FastqRead(&fq)) {
            for (b=fq.read.qual; *b != '\n'; ++b) {
                min = *b < min ? *b : min;
                max = *b > max ? *b : max;
            }
        }
        else break;
    } gzclose(fq.in);

    if (min < 64) return Phrd33;
    if (max > 75) return Phrd64;
}

/*! @funciton: read length check
 *   @parmeters:
 *   @    infile     the pointer to the input fastq file [char *]
 *   @return:
 *   @    maxlen     the maximum read length
*/
int ReadLenCheck(char *infile)
{
    int maxlen=0;
    fastq_t fq;

    err_gzopen(fq.in, infile, "r");
    for (int i=0; i < 1000; ++i) {
        if (FastqRead(&fq)) {
            int curlen = strlen(fq.read.seq) - 1;
            maxlen = curlen > maxlen ? curlen : maxlen;
        }
        else break;
    } gzclose(fq.in);

    if (maxlen >= FQLINE) {
        fprintf(stderr, "[Error] Read length is longer than defined maximum read length(%d)\n", FQLINE);
        exit(-1);
    }
      
    return maxlen;
}

/*! @funciton: calculate the average quality
 *   @parmeters:
 *   @    qual       the pointer to the quality [char *]
 *   @    phred      the Phred encoding [int]
 *   @return:
 *   @    mean quality
*/
float MeanQuality(char *qual, int phred)
{
    int q_sum=0, num=0;
    int ascii[2] = {33, 64};

    for (qual; *qual; ++qual) {
        if (*qual != '\r' && *qual != '\n') {
            q_sum += *qual - ascii[phred];
            ++num;
        }
    }
    
    return (float)q_sum/num;
}


/*! @funciton: read a fastq group each time [seqname, seq, +, qual]
 *   @parmeters:
 *   @    fq         the pointer to the fastq_t structure [fastq_t *]
 *   @return:
 *   @    void
*/
int FastqRead(fastq_t *fq)
{
    char buf[BUFLINE];

    for (int i=0; i < 4; i++) {
        if (gzgets(fq->in, buf, BUFLINE)) {
            switch (i) {
                case 0: strcpy(fq->read.name, buf); break;
                case 1: strcpy(fq->read.seq,  buf); break;
                case 2: strcpy(fq->read.mark, "+\n"); break;
                case 3: strcpy(fq->read.qual, buf);
            }
        }
        else return 0;
    }
    return 1;
}


/*! @funciton: write the fastq group from cache
 *   @parmeters:
 *   @    fq         the pointer to the fastq_t structure [fastq_t *]
 *   @return:
 *   @    void
*/
void FastqWrite(fastq_t *fq)
{
    read_t *read;
   
    read = fq->cache;
    for (int i=0; i < fq->bufnum; i++) {
        for (int j=0; j < 4; j++) {
            switch (j) {
                case 0: fputs(read[i].name, fq->out); break;
                case 1: fputs(read[i].seq,  fq->out); break;
                case 2: fputs(read[i].mark, fq->out); break;
                case 3: fputs(read[i].qual, fq->out);
            }
        }
    } fq->bufnum = 0;
    return ;
}

