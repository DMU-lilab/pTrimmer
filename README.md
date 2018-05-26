pTrimmer
=========================
Used to trim off the primer sequence from amplicon fastq file


__PROGRAM: pTrimmer__<br>
__VERSION: 1.3.1__<br>
__PLATFORM: Linux and Windows__<br>
__COMPILER: gcc-4.8.5__<br>
__AUTHOR: xiaolong zhang__<br>
__EMAIL: xiaolongzhang2015@163.com__<br>
__DATE:   2017-09-21__<br>
__UDDATE: 2018-05-09__<br>
__DEPENDENCE__<br>
* zlib-1.2.7<br>
#### NOTE
* The first you need to do are confirming the libraries above have been installed.<br />
* And the gcc compiler should be available on your server or laptop.<br /><br />
* The program could run on a standard dual core laptops with 8 GB of RAM on both windows(win7 or win10) and linux(centos or ubuntu).<br /><br />


Description
=========================
* The program is used to trim off the primer sequence of the target sequencing at both 5'(forward primer) and 3'(reverse complement primer). <br>
* Both k-mer indexing alogrithm and dynamic programing alogrithm were performed to trim off the primer sequencemakes.<br>
* The performing of k-mer (seed and extend) alogrithm makes it possible to deal with __thousands of amplicon primer pairs__ at the same time.<br>
* Compared with other kinds of tools, this program could trim the primer sequence off directly from the fastq file, which could save you a lot of time.<br>
* There only have 250 reads in the example fastq file, which result in a higher mismatch ratio. But in the general amplicon data, the program will has a good performance and lower mismatch ratio.


Building
=========================

See INSTALL for complete details.


Usage
========================
      Options:
         -h|--help        print help infomation
         -l|--keep        keep the complete reads if can't locate primer
                          sequence [default: discard the reads]
         -s|--seqtype     [required] the sequencing type [single|pair]
         -a|--ampfile     [required] input amplicon file [.txt]
         -f|--read1       [required] read1(forward) for fastq file [.fq|.gz]
         -r|--read2       [optional] read2(reverse) for paired-end seqtype [.fq|.gz]
         -o|--outdir      [required] output directory for trimed fastq file [dir]
         -q|--minqual     [optional] the minimum average quality to keep after triming [20]
         -k|--kmer        [optional] the kmer lenght for indexing [8]
         -m|--mismatch    [optional] the maxmum mismatch for primer seq [3]


Option
========================
#### \[-h|--help]
      Print out the help infomation

#### \[-l|--keep]
      Whether keep the original complete reads if failed to locate the 
      primer sequence [default: discard the reads]

#### \[-s|--seqtype]
      The sequencing type of input fastq file, which could be single-end 
      or paired-end [single or pair]

#### \[-a|--ampfile]
      The path to amplicon file with text format.
      There is an example in the "Example" directory, which will give you a guidance
      about how to specify each field of the amplicon file.
      [FIELDS]
            1. forwardprimer: the forward primer sequecne (5') [required]
            2. reverseprimer: the reverse primer sequence (5') [required]
            3. insertlength: the insert length between the primer pair [required]
            4. AuxInfo: the auxiliary information that used to describe 
               the primer pair [optional]
      eg. data_amplicon.txt

#### \[-f|--read1]
      The read1 of the fastq file or gziped fastq file [must have an marker of "_R"]
      eg. data_R1.fq or data_R1.fq.gz

#### \[-r|--read2]
      The read2 of the fastq file or gziped fastq file [must have an marker of "_R"]
      eg. data_R2.fq or data_R2.fq.gz

#### \[-o|--outdir]
      The output directory of the trimed fastq file

#### \[-q|--minqual]
      The minimum average quality to keep after triming. The program will automatically
      detect the quality coding(Phrd+33 or Phrd+64) [default: 20]

#### \[-k|--kmer]
      The kmer length for primer indexing, a larger kmer will perform an accurate locating
      to the primer sequence, which also means a lower posibility to find the primer 
      sequence. [8 is recommended]

#### \[-m|--mismatch]
      the maxmum mismatch allowed in the  primer matching [ mismatch <= 3 is recommended]


OutPut
=========================
The program will generate 3 files for paired-end:
1. \*_trim_R1.fq
2. \*_trim_R2.fq
3. Summary.ampcount

or 2 files for single-end:
1. \*_trim_R1.fq
2. Summary.ampcount

Description:
1. \*_trim_R1.fq and \*_trim_R2.fq are the fastq file after primer trimed!<br>
2. The field of "AmpCount" in the Summary.ampcount file indicate the reads belong to the amplicon.<br>
