/* function: parse the comand line parmeters */

#include "fastq.h"
#include "utils.h"

void Usage(void)
{
    char *usage =
        "\nUsage: pTrimmer [options]\n"
        "Version: 1.3.1\n"
        "\n"
        "Options:\n"
        "       -h|--help        print help infomation\n"
        "       -l|--keep        keep the complete reads if can't locate primer\n" 
        "                        sequence [default: discard the reads]\n"
        "       -s|--seqtype     [required] the sequencing type [single|pair]\n"
        "       -a|--ampfile     [required] input amplicon file [.txt]\n"
        "       -f|--read1       [required] read1(forward) for fastq file [.fq|.gz]\n"
        "       -r|--read2       [optional] read2(reverse) for paired-end seqtype [.fq|.gz]\n"
        "       -o|--outdir      [required] output directory for trimed fastq file [dir]\n"
        "       -q|--minqual     [optional] the minimum average quality to keep after triming [20]\n"
        "       -k|--kmer        [optional] the kmer lenght for indexing [8]\n"
        "       -m|--mismatch    [optional] the maxmum mismatch for primer seq [3]\n\n";

    fprintf(stderr, "%s", usage); 
    exit(-1);
}

static const struct option long_options[] =
{
    { "help", no_argument, NULL, 'h' },
    { "keep", no_argument, NULL, 'l' },
    { "seqtype", required_argument, NULL, 's' },
    { "ampfile", required_argument, NULL, 'a' },
    { "read1", required_argument, NULL, 'f' },
    { "read2", required_argument, NULL, 'r' },
    { "outdir", required_argument, NULL, 'o' },
    { "minqual", required_argument, NULL, 'q' },
    { "kmer", required_argument, NULL, 'k' },
    { "mismatch", required_argument, NULL, 'm' },
    { NULL, 0, NULL, 0 }
};

static void ArgInit(arg_t *Arg)
{
    Arg->keep = 0;
    Arg->seqtype = -1;
    Arg->minqual = 20;
    Arg->kmer = 8;
    Arg->mismatch = 3;
}

arg_t *ParseOpt( int argc, char **argv )
{
    int opt =0, opterr =0;
    arg_t *Arg;
    
    err_calloc(Arg, 1, arg_t);
    ArgInit(Arg);
    while ( (opt = getopt_long(argc, argv, "s:a:f:r:o:q:k:m:hl", long_options, NULL)) != -1 )
    {
        switch (opt) {
            case 'h': Arg->help = 1; break;
            case 'l': Arg->keep = 1; break;
            case 's': if (!strcmp(optarg, "single")) Arg->seqtype = SE;
                      else if (!strcmp(optarg, "pair")) Arg->seqtype = PE;
                      else Arg->seqtype = -1; break;
            case 'a': strcpy(Arg->ampfile, optarg); break;
            case 'f': strcpy(Arg->read1, optarg); break;
            case 'r': strcpy(Arg->read2, optarg); break;
            case 'o': strcpy(Arg->outdir, optarg); break;
            case 'q': Arg->minqual = atoi(optarg); break;
            case 'k': Arg->kmer = atoi(optarg); break;
            case 'm': Arg->mismatch = atoi(optarg); break;
            case '?': fprintf(stderr, \
                            "[Err::%s::%d]  Option error occour!.\n", __func__, __LINE__);
                      Arg->help = 1;
        }
    } 
    if (!Arg->ampfile[0] || !Arg->read1[0] || (Arg->seqtype == -1) || !Arg->outdir[0]) {
        fprintf(stderr, \
            "[Err::%s::%d]  Please give the [requied] parmeters!\n", __func__, __LINE__);
        Arg->help = 1;
    }

    return Arg;
}

