#!/usr/bin/env python
# -*- coding: utf-8 -*-

#*************************************************************************
#    > File Name: Bed2Amplicon.py
#    > Author: xlzh
#    > Mail: xiaolongzhang2015@163.com 
#    > Created Time: 2021年06月08日 星期二 10时12分15秒
#*************************************************************************

import sys
import os.path


class Fasta(object):
    def __init__(self, fa_file):
        self.fa_file = fa_file
        self.__fa_fp = self.__fa_init()
        self.__fai_dict = self.__fai_load()

    def __fa_init(self):
        ''' Func: open the fasta handle
        '''
        if not os.path.exists(self.fa_file):
            sys.stderr.write('Err: No such file: %s\n' %(self.fa_file))
            sys.exit(-1)

        return open(self.fa_file, 'r')

    def __fai_load(self):
        ''' Func: load the index of fasta file
        '''
        if not os.path.exists(self.fa_file + '.fai'):
            sys.stderr.write('Err: No such index file: %s\n' %(self.fa_file+'.fai'))
            sys.exit(-1)

        fai_dict = {} # {'chr1':(offset, b_len, l_len), 'chr2':(), ...}

        fai_fp = open(self.fa_file+'.fai', 'r')
        for line in fai_fp:
            l = line.split()
            fai_dict[l[0]] = (int(l[2]), int(l[3]), int(l[4]))

        return fai_dict

    def __read_base(self):
        ''' Func: read 1-base each time
        '''
        base = self.__fa_fp.read(1)
        if base == '': return 0

        return base

    def get_seq(self, chrom, start, end):
        ''' Func: get the fasta sequence from start to end
        '''
        seq = []

        if (chrom not in self.__fai_dict) or (end < start):
            sys.stderr.write('Warning: Invaild input: %s:%d-%d\n' %(chrom,start,end))
            return ''

        idx = self.__fai_dict[chrom]
        d_offset = idx[0] + (start-1) + (start-1)/idx[1]*(idx[2]-idx[1])

        self.__fa_fp.seek(d_offset) # put the file handle to the start position
        ch, seq_len = 1, end-start+1

        while (ch and seq_len):
            ch = self.__read_base()
            if ch not in ['\r', '\n', 0]: seq.append(ch); seq_len -= 1

        return ''.join(seq)


def get_reverse_complement(primer_seq):
    ''' Func: get the reverse complement sequence
    '''
    new_seq_list = []
    Table = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

    for s in primer_seq[::-1]:
        new_seq_list.append(Table[s])

    return ''.join(new_seq_list)


def show_error_message():
    sys.stderr.write("Error: invalid bed file format!\n")
    sys.stderr.write("bed file must be have 4 column and the 4-th column must has a marker of 'R' or 'F'\n")
    sys.stderr.write("to indicate reverse and forward primer\n")
    sys.stderr.write("eg. NC_045512.2\t33\t56\tAmplicon1_F\n")
    sys.stderr.write("    NC_045512.2\t88\t102\tAmplicon1_R\n")
    sys.exit(-1)


def read_bed_file(bed_file, fasta_obj):
    ''' Func: read input bed file
        bed_dict = {
            'amplicon1': [(forward record), (reverse record)], 
            'amplicon2': [(chr1, 3, 26), (chr1, 88, 102)]
            ...
        }
        Note: the last field must has a marker of 'R' or 'F' to indicate reverst and forward strand
              position is 0-based
    '''
    bed_dict = {}
    bed_fp = open(bed_file, "r")

    for line in bed_fp:
        llist = line.rstrip().split()
        if len(llist) < 4: show_error_message()

        amplicon = llist[3][:-1]
        strand_marker = llist[3][-1]
        if strand_marker not in ['F', 'R']: show_error_message()

        primer_seq = fasta_obj.get_seq(llist[0], int(llist[1])+1, int(llist[2])+1)
        if amplicon not in bed_dict: bed_dict[amplicon] = [None, None]

        if strand_marker == 'F':  # forward strand
            bed_dict[amplicon][0] = (primer_seq, int(llist[2]))
        else:  # reverse strand and has a marker of 'R'
            bed_dict[amplicon][1] = (get_reverse_complement(primer_seq), int(llist[1]))

    return bed_dict


def write_amplicon_file(amplicon_file, bed_dict):
    ''' Func: write the amplicon file out
    '''
    amplicon_fp = open(amplicon_file, 'w')
    amplicon_fp.write("#ForwardPrimer\tReversePrimer\tInsertLength\tAdditionInfo\n")

    for amplicon in bed_dict:
        bed_tuple = bed_dict[amplicon]
        if not bed_tuple[0] or not bed_tuple[1]:  
            continue # didn't find mathed bed record with different strand marker 'R' or 'F'

        insert_len = bed_tuple[1][1] - bed_tuple[0][1] + 1
        amplicon_fp.write("%s\t%s\t%d\t%s\n" % (bed_tuple[0][0], bed_tuple[1][0], insert_len, amplicon))

    amplicon_fp.close()


def main():
    args = sys.argv

    if len(args) != 4:
        sys.stderr.write("usage: python Bed2Amplicon.py <genome.fa> <input.bed> <out_amplicon.txt>\n")
        sys.exit(-1)

    refer_file = args[1]  # reference genome (.fasta)
    bed_file = args[2]   # bed file with 0-based reference position
    amplicon_file = args[3]  # converted amplicon primer file

    # refer_file = "/Users/xlzh/Downloads/PTrimmer/Data/NC_045512.2.fas"
    # bed_file = "/Users/xlzh/Downloads/PTrimmer/Data/sarscov2_v2_primers_swift.bed"
    # out_file = "/Users/xlzh/Downloads/PTrimmer/Data/amplicon.txt"

    fasta_obj = Fasta(refer_file)
    bed_dict = read_bed_file(bed_file, fasta_obj)
    write_amplicon_file(amplicon_file, bed_dict)


if __name__ == '__main__':
    main()

