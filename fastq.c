#include "fastq.h"
#include "utils.h"


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
    fq->size = 0;
    fq->max = BUFF_NUM;
    memset(&fq->read, 0, sizeof(read_t));

    if (type == READ1) {
        fq->in_hd = gz_stream_open(args->read1, "r");
        fq->out_hd = gz_stream_open(args->trim1, "w");
    }
    else {
        fq->in_hd = gz_stream_open(args->read2, "r");
        fq->out_hd = gz_stream_open(args->trim2, "w");
    }

    /* allocate memory for the read cache */
    err_calloc(fq->cache, fq->max, read_t);
}


/*! @funciton: Phred encoding check
 *   @parmeters:
 *   @    infile     the pointer to the input fastq file [char *]
 *   @return:
 *   @    Phred33[0] or Phred64[1]
*/
int PhredCheck(char *infile)
{
    int min=128, max=0;
    read_t read;

    GzStream *gz_hd = gz_stream_open(infile, "r");
    for (int i=0; i < 1000; i++) {
        if (FastqGetRead(gz_hd, &read) == 0)
            break;  /* EOF of the fastq file */

        kstring_t *qual = &read.qual;
        int qual_len = (int)qual->l - 2;  /* skip the '\r' and '\n' */

        for (int j=0; j < qual_len; j++) {
            if (qual->s[j] < min) min = (int)qual->s[j];
            if (qual->s[j] > max) max = (int)qual->s[j];
        }
    }
    gz_stream_destroy(gz_hd);
    FastqFreeRead(&read);

    if (min < 64) return Phrd33;
    if (max > 75) return Phrd64;

    return Phrd33;
}

/*! @funciton: read length check
 *   @parmeters:
 *   @    infile     the pointer to the input fastq file [char *]
 *   @return:
 *   @    max_len    the maximum read length
*/
int ReadLenCheck(char *infile)
{
    int max_len = 0;
    read_t read;

    GzStream *gz_hd = gz_stream_open(infile, "r");
    for (int i=0; i < 1000; i++) {
        if (FastqGetRead(gz_hd, &read) == 0)
            break;  /* EOF of the fastq file */

        int seq_len = (int)read.seq.l;
        if (seq_len > max_len) max_len = seq_len;
    }
    gz_stream_destroy(gz_hd);
    FastqFreeRead(&read);
      
    return max_len;
}

/*! @funciton: calculate the average quality
 *   @parmeters:
 *   @    qual       the pointer to the quality [char *]
 *   @    phred      the Phred encoding [int]
 *   @return:
 *   @    mean quality
*/
float MeanQuality(kstring_t *qual_obj, int phred)
{
    int q_sum=0, num=0;
    int ascii[2] = {33, 64};

    for (char *qual=qual_obj->s; *qual; ++qual) {
        if (*qual != '\r' && *qual != '\n') {
            q_sum += *qual - ascii[phred];
            ++num;
        }
    }

    return (float)q_sum / (float)num;
}


static int FastqReadCore(GzStream *gz, read_t *read)
{
    int ret;

    ret = gz_read_util(gz, '\n', &read->name);  /* get the read name*/
    switch (ret) {
        case -1: return -1;  /* unexpected end of the file */
        case  0: return 0;  /* end of the file or empty file */
        default: break;  /* normal reading the file (ret==1) */
    }

    ret = gz_read_util(gz, '\n', &read->seq);  /* get the read sequence */
    switch (ret) {
        case -1: return -1;  /* unexpected end of the file */
        case  0: return -2;  /* incomplete fastq reads, since only one line is readed */
        default: break;  /* normal reading the file (ret==1) */
    }

    ret = gz_read_util(gz, '\n', &read->comment); /* get the read comment, usually is '+' */
    switch (ret) {
        case -1: return -1;  /* unexpected end of the file */
        case  0: return -2;  /* incomplete fastq reads, since only two lines are readed */
        default: break;  /* normal reading the file (ret==1) */
    }

    ret = gz_read_util(gz, '\n', &read->qual); /* get the read quality */
    switch (ret) {
        case -1: return -1;  /* unexpected end of the file */
        case  0: return -2;  /* incomplete fastq reads, since only three lines are readed */
        default: break;  /* normal reading the file (ret==1) */
    }

    return ret;  /* means read status is OK */
}


/*! @funciton: read the input fastq reads
 *   @parmeters:
 *   @    gz_hd         the pointer to the input stream
 *   @    read          the pointer to the dest read
 *   @return:
 *   @    status        0: EOF; 1: OK
*/
int FastqGetRead(GzStream *gz_hd, read_t *read)
{
    int ret;

    ret = FastqReadCore(gz_hd, read);
    switch (ret) {
        case 0:
            return 0;  /* the really EOF */
        case -1:
            fprintf(stderr, "[SysError:fastq_cache_read] unexpected end (truncated) of fastq file is detected!\n");
            exit(-1);
        case -2:
            fprintf(stderr, "[SysErr:fastq_cache_read] incomplete fastq read is detected!\n");
            exit(-1);
        default:
            break; /* normal reading of the fastq read (ret==1) */
    }

    return 1;
}


/*! @funciton: write the trimmed fastq from cache
 *   @parmeters:
 *   @    fq         the pointer to the fastq_t structure [fastq_t *]
 *   @return:
 *   @    void
*/
void FastqWriteCache(fastq_t *fq)
{
    read_t *read;

    for (int i=0; i < fq->size; i++) {
        read = &fq->cache[i];

        gz_write_util(fq->out_hd, &read->name);
        gz_write_util(fq->out_hd, &read->seq);
        gz_write_util(fq->out_hd, &read->comment);
        gz_write_util(fq->out_hd, &read->qual);
    }

    fq->size = 0;
}


void FastqDestroy(fastq_t *fq)
{
    gz_stream_destroy(fq->in_hd);
    gz_stream_destroy(fq->out_hd);
    if (fq->cache) free(fq->cache);
}


void FastqFreeRead(read_t *read)
{
    /* destroy read name */
    if (read->name.s) {
        free(read->name.s);
        read->name.s = NULL; read->name.l = read->name.m = 0;
    }

    /* destroy read seq */
    if (read->seq.s) {
        free(read->seq.s);
        read->seq.s = NULL; read->seq.l = read->seq.m = 0;
    }

    /* destroy read comment */
    if (read->comment.s) {
        free(read->comment.s);
        read->comment.s = NULL; read->comment.l = read->comment.m = 0;
    }

    /* destroy read quality */
    if (read->qual.s) {
        free(read->qual.s);
        read->qual.s = NULL;
        read->qual.l = read->qual.m = 0;
    }
}
