/* function: parse the comand line parameters */

#include "parse.h"
#include "utils.h"
#include "version.h"


void Usage(void)
{
    char *usage =
        "\nUsage: pTrimmer [options]\n"
        "\n"
        "Options:\n"
        "       -h|--help        print help information\n"
        "       -l|--keep        keep the complete reads if failed to locate primer\n" 
        "                        sequence [default: discard the reads]\n"
        "       -t|--seqtype     [required] the sequencing type [single|pair]\n"
        "       -a|--ampfile     [required] input amplicon file [.txt]\n"
        "       -f|--read1       [required] read1(forward) for fastq file [.fq|.gz]\n"
        "       -d|--trim1       [required] the trimmed read1 of fastq file\n"
        "       -r|--read2       [optional] read2(reverse) for fastq file (paired-end seqtype) [.fq|.gz]\n"
        "       -e|--trim2       [optional] the trimmed read2 of fastq file (paired-end seqtype)\n"
        "       -z|--gzip        [optional] output trimmed fastq file in Gzip format\n"
        "       -i|--info        [optional] add the primer information for each trimmed read\n"
        "       -s|--summary     [optional] the trimming information of each amplicon [default: Summary.ampcount]\n"
        "       -q|--minqual     [optional] the minimum average quality to keep after trimming [20]\n"
        "       -k|--kmer        [optional] the kmer length for indexing [8]\n"
        "       -m|--mismatch    [optional] the maximum mismatch for primer seq [3]\n\n";

    fprintf(stderr, "Program: pTrimmer (v%s)\n", PTRIMMER_VERSION_STRING);
    fprintf(stderr, "CreateDate: %s\n", PTRIMMER_CREATE_DATE);
    fprintf(stderr, "UpdateDate: %s\n", PTRIMMER_UPDATE_DATE);
    fprintf(stderr, "Author: XiaolongZhang (xiaolongzhang2015@163.com)\n");
    fprintf(stderr, "%s", usage); 
    exit(-1);
}

static const struct option long_options[] =
{
    { "help", no_argument, NULL, 'h' },
    { "keep", no_argument, NULL, 'l' },
    { "seqtype", required_argument, NULL, 't' },
    { "ampfile", required_argument, NULL, 'a' },
    { "read1", required_argument, NULL, 'f' },
    { "trim1", required_argument, NULL, 'd' },
    { "read2", required_argument, NULL, 'r' },
    { "trim2", required_argument, NULL, 'e' },
    { "gzip", no_argument, NULL, 'z' },
    { "info", no_argument, NULL, 'i' },
    { "summary", required_argument, NULL, 's' },
    { "minqual", required_argument, NULL, 'q' },
    { "kmer", required_argument, NULL, 'k' },
    { "mismatch", required_argument, NULL, 'm' },
    { NULL, 0, NULL, 0 }
};

static void ArgInit(arg_t *Arg)
{
    Arg->keep = 0;
    Arg->gzip = 0;
    Arg->info = 0;
    Arg->seqtype = -1;
    Arg->minqual = 20;
    Arg->kmer = 8;
    Arg->mismatch = 3;
}

arg_t *ParseOpt( int argc, char **argv )
{
    int opt;
    arg_t *Arg;
    
    err_calloc(Arg, 1, arg_t);
    ArgInit(Arg);
    while ( (opt = getopt_long(argc, argv, "t:a:f:d:r:e:s:q:k:m:hlzi", long_options, NULL)) != -1 )
    {
        switch (opt) {
            case 'h': Arg->help = 1; break;
            case 'l': Arg->keep = 1; break;
            case 'z': Arg->gzip = 1; break;
            case 'i': Arg->info = 1; break;
            case 't': if (!strcmp(optarg, "single")) Arg->seqtype = SE;
                      else if (!strcmp(optarg, "pair")) Arg->seqtype = PE;
                      else Arg->seqtype = -1; break;
            case 'a': strcpy(Arg->ampfile, optarg); break;
            case 'f': strcpy(Arg->read1, optarg); break;
            case 'd': strcpy(Arg->trim1, optarg); break;
            case 'r': strcpy(Arg->read2, optarg); break;
            case 'e': strcpy(Arg->trim2, optarg); break;
            case 's': strcpy(Arg->summary, optarg); break;
            case 'q': Arg->minqual = atoi(optarg); break;
            case 'k': Arg->kmer = atoi(optarg); break;
            case 'm': Arg->mismatch = atoi(optarg); break;
            default:
                fprintf(stderr, "[Err::%s::%d] Option error occurred!.\n", __func__, __LINE__);
                Arg->help = 1;
        }
    } 
    if (!Arg->ampfile[0] || !Arg->read1[0] || (Arg->seqtype == -1) || !Arg->trim1[0]) {
        fprintf(stderr, "[Err::%s::%d] Please give the [required] parameters!\n", __func__, __LINE__);
        Arg->help = 1;
    }
    /* check the suffix of the trimmed fastq file */
    if (Arg->gzip && !strstr(Arg->trim1, ".gz")) {
        strcat(Arg->trim1, ".gz");  /* add suffix of .gz to the file */
    }

    if (Arg->seqtype == PE) {
        if (!Arg->read2[0] || !Arg->trim2[0]) {
            fprintf(stderr, "[Err::%s::%d] Please give the parameters (read2 and trim2) in pair-end mode!\n", __func__, __LINE__);
            Arg->help = 1;
        }
        /* check the suffix of the trimmed fastq file */
        if (Arg->gzip && !strstr(Arg->trim2, ".gz")) {
            strcat(Arg->trim2, ".gz");  /* add suffix of .gz to the file */
        }
    }

    return Arg;
}

