/*************************************************************************
    > File Name: fileio.h
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com
    > Created Time: 2024年01月10日 星期三 12时30分15秒
 ************************************************************************/

#include "fileio.h"


/* check whether the str is ends with sub.
 * e.g. get_str_ends('/home/data/test.fq.gz', '.gz') -> true
 * */
static int get_str_ends(char *str, char *sub)
{
    int i = (int)strlen(str) - 1;
    int j = (int)strlen(sub) - 1;

    if(i < j) return 0;

    for(; i>=0 && j>=0; i--,j--){
        if(str[i] != sub[j]) return 0;
    }
    return 1;
}


/* get the basename of the given file path */
static char *get_path_basename(char *file_path)
{
#ifdef _WIN32
    #define PATH_DELIM 92  /* '\' */
#else
    #define PATH_DELIM 47  /* '/' */
#endif

    char *start;

    start = strrchr(file_path, PATH_DELIM);
    if (start != NULL) return start + 1;

    return file_path;
}


GzStream *gz_stream_open(char *file, char *mode)
{
    if(strcmp(mode, "w") != 0 && strcmp(mode, "r") != 0) {
        fprintf(stderr, "[SysError:gz_stream_open:01]: operate mode(%s) error, it should be \"w\" or \"r\".\n", mode);
        exit(0);
    }

    size_t n_items;
    GzStream *gz = calloc(1, sizeof(GzStream));
    char *err_fn = get_path_basename(file);
    if(mode[0]=='r') {
        if(get_str_ends(file, ".gz")) {
            unsigned char buf[4];
            FILE *f = fopen(file, "rb");
            if (f==NULL) {
                fprintf(stderr, "[SysError:gz_stream_open:02]: can not open gzip file of (%s)!\n", err_fn);
                exit(0);
            }
            n_items = fread(buf, 1, 4, f);  /* read the magic number of the gzip file */
            if (n_items != 4) {
                fprintf(stderr, "[SysError:gz_stream_open:03]: truncated file header detected in file of (%s)!\n", err_fn);
                exit(0);
            }
            if (buf[0]==0x1f && buf[1]==0x8b) {
                gz->gz_fp=gzopen(file, "r");
                if(gz->gz_fp <= 0){
                    fprintf(stderr, "[SysError:gz_stream_open:04]: failed to open file of (%s)!\n", err_fn);
                    exit(0);
                }
            }
            else {
                fprintf(stderr, "[SysError:gz_stream_open:04]: the file (%s) is not a gzip file!\n", err_fn);
                exit(0);
            }
            fclose(f);
        }
        else {  /* normal (uncompressed) fastq file */
            unsigned char buf[4];
            FILE *f=fopen(file, "rb");
            if(f==NULL) {
                fprintf(stderr, "[SysError:gz_stream_open:08]: can not open fastq file of (%s)!\n", err_fn);
                exit(0);
            }
            n_items = fread(buf, 1, 4, f);
            if (n_items != 4) {
                fprintf(stderr, "[SysError:gz_stream_open:09]: truncated file header detected in file of (%s)!\n", err_fn);
                exit(0);
            }
            if (buf[0] != '@') {
                fprintf(stderr, "[SysError:gz_stream_open:10]: the file of (%s) is neither a real gzip(.gz)/bzip2(.bz2) file nor a stander fastq file!\n", err_fn);
                exit(0);
            }
            gz->gz_fp=gzopen(file, "r");
            fclose(f);
        }

        gz->is_write=0;
    }
    else {  /* open a gz or normal file for writing */
        if(get_str_ends(file, ".gz")) {
            gz->gz_fp=gzopen(file, "w");
            if(gz->gz_fp<=0){
                fprintf(stderr, "[SysError:gz_stream_open:11]: failed to create file of (%s) in .gz format!\n", err_fn);
                exit(0);
            }
        }
        else{
            gz->out_fp=fopen(file, "w");
            if(gz->out_fp<=0){
                fprintf(stderr, "[SysError:gz_stream_open:14]: failed to create file of (%s)!\n", err_fn);
                exit(0);
            }
        }
        gz->is_write=1;
    }
    gz->buf = calloc(GZ_BUFF_SIZE, sizeof(char));

    return gz;
}


int gz_read_util(GzStream *gz, char delimiter, kstring_t *ks_str)
{
    int len = 0;
    char c;

    if(gz->is_eof && gz->begin>=gz->end) return 0;  /* end of the file */

    do {
        if (gz->begin >= gz->end){ /* gz->buf is full or the first time to read */
            gz->begin = 0;
            if (gz->gz_fp) {
                gz->end = gzread(gz->gz_fp, gz->buf, GZ_BUFF_SIZE);
                gzerror(gz->gz_fp, &(gz->bzerror));
                if (gz->bzerror < 0) return -1; /* truncated file detected */
            }
            if (gz->end < GZ_BUFF_SIZE) gz->is_eof = 1;
        }
        while (gz->begin < gz->end) { /* copy the str to the user-provided ks_str */
            c = gz->buf[gz->begin++];
            if (len+1 >= ks_str->m) {  /* specially for space of delimiter */
                ks_str->m = ks_str->m ? ks_str->m<<1 : 256;
                if (ks_str->m > READ_LEN_MAX) {
                    fprintf(stderr, "[SysError:gz_read_util:014] read length can not be longer than %d!\n", READ_LEN_MAX);
                    exit(-1);
                }
                ks_str->s = (char *)realloc(ks_str->s, ks_str->m * sizeof(char));
                if (!ks_str->s) {
                    fprintf(stderr, "[SysError:gz_read_util:015] failed to reallocated memory!\n");
                    exit(-1);
                }
            }
            if (c == delimiter) {  /* add the delimiter('\n') to the end of the string */
                ks_str->s[len++] = delimiter; ks_str->s[len] = '\0';
                ks_str->l = len;
                return 1;
            }
            ks_str->s[len++] = c;
        }
    } while (!gz->is_eof);

    return 0;  /* end of the file (EOF) */
}


void gz_write_util(GzStream *gz, kstring_t *ks_str)
{
    char *str = ks_str->s;
    size_t str_len = ks_str->l;

    while(1) {
        int last = GZ_BUFF_SIZE - gz->begin;

        if(str_len <= last) {
            memcpy(gz->buf+gz->begin, str, str_len);
            gz->begin += (int)str_len;
            return;
        }
        memcpy(gz->buf+gz->begin, str, last);
        if(gz->gz_fp) gzwrite(gz->gz_fp, gz->buf, GZ_BUFF_SIZE);
        else fwrite(gz->buf, 1, GZ_BUFF_SIZE, gz->out_fp);
        
        gz->begin = 0;
        str += last;
        str_len -= last;
    }
}


void gz_stream_destroy(GzStream *gz)
{
    if (gz->is_write){  /* write gz or gz2 file */
        if(gz->gz_fp){
            gzwrite(gz->gz_fp, gz->buf, gz->begin);
            gzclose(gz->gz_fp);
        }
        else {  /* normal file writing */
            fwrite(gz->buf, 1, gz->begin, gz->out_fp);
            fclose(gz->out_fp);
        }
    }
    else {  /* read gz or bz2 file */
        if (gz->gz_fp) gzclose(gz->gz_fp);
    }
    if(gz->buf) free(gz->buf);
    free(gz);
}
