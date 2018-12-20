#!/usr/bin/python

import pdb
import sys
from Bio import SeqIO
import gzip
import random


class CreateRead(object):
    def __init__(self):
        self.qname = ''
        self.seq = ''
        self.black = '+'
        self.quality = ''

    def qnameassign(self,qname):
        self.qname = qname

    def seqassign(self,seq):
        self.seq = seq

    def qualityassign(self,quality):
        self.quality = quality

class FqWrite(object):
    def __init__(self, outBase):
        self.file1 = gzip.open(outBase + "_R1.fq.gz", "a")
        self.file2 = gzip.open(outBase + "_R2.fq.gz", "a")

    def write(self, read1, read2):
        self.file1.write("%s\n%s\n" %(read1.qname,read1.seq))
        self.file1.write("%s\n%s\n" %(read1.black,read1.quality))
        self.file2.write("%s\n%s\n" %(read2.qname,read2.seq))
        self.file2.write("%s\n%s\n" %(read2.black,read2.quality))

    def close(self):
        self.file1.close()
        self.file2.close()
        
def rev_complement(seq):
        seq = seq.upper()
        basecomplement = {'A':'T',
                          'C':'G',
                          'G':'C',
                          'T':'A',
                          '-':'-',
                          'N':'N'}
        letters = list(seq)
        letters = [basecomplement[base] for base in letters]
        complement = (''.join(letters)) #gives the complement of the bases in list letters
        return complement[::-1]

def GetGenome(inputfile):
    genome_dict = dict()
    fasta_seq = SeqIO.parse(inputfile, "fasta")

    for fasta in fasta_seq:
        key = fasta.id
        value = fasta.seq
        if key not in genome_dict:
            genome_dict[key] = value
        else:
            sys.stderr.write("Err: duplicate chromosome occured in your genome file!!!\n")
            sys.exit(-1)

    return genome_dict

def noise(rate, primseq):
    ''' M -> mismatch; I -> insert; D -> deletion
    '''
    err_type, n = ['M','I','D'], random.uniform(0,1)
    if n > rate:
        return (0, primseq)

    plen = len(primseq)
    idx, n1 = random.randint(0, plen-1), random.randint(0, 2)

    if err_type[n1] == 'M':
        seq = primseq[:idx] + 'N' + primseq[idx+1:plen]
        return (0, seq)

    if err_type[n1] == 'I':
        seq = primseq[:idx] + 'N' + primseq[idx:]
        return (1, seq)

    if err_type[n1] == 'D':
        seq = primseq[:idx] + primseq[idx+1:]
        return (-1, seq)

def main():
    Usage = [ "Function: simulate multiple amplicon sequencing data\n\n",
              "Usage: python MASimulator.py <ref.fa> <primer.csv> <out_prefix> <in-target.depth> <off-target.depth> <noise_ratio>\n\n",
              "Note:\n",
              "   1. ref.fa: reference genome, eg. hg19.fa\n",
              "   2. primer.csv: contain chromosome and position of amplicon, eg. CfDNA_Full.csv\n",
              "   3. target.depth: depth of in-target NGS reads\n",
              "   4. off-target.depth: depth of off-target NGS reads\n",
              "   5. noise_ratio: introduced mismatch/insertation/deletion ratio in the primer regions\n"
            ]
             

    args = sys.argv

    if len(args) == 1:
        sys.stderr.write(''.join(Usage))
        sys.exit(-1)

    genome_dict = GetGenome(args[1])
    ampliconcsv = open(args[2], "r")
    Fp = FqWrite(args[3])

    intargat_reads_depth, oftarget_reads_depth = int(args[4]), int(args[5])
    noise_rate, pos, amplicon_dict = float(args[6]), 0, dict()

    for line in ampliconcsv.readlines():
        if line.startswith('#'):
            continue

        lines = line.rstrip().split(",")
        pos += 1
        amplicon_chrom = lines[7]
        amplicon_start, amplicon_end = int(lines[8]), int(lines[11])
        insert_start, insert_end = int(lines[9]), int(lines[10])

        insert_len = insert_end - insert_start 
        forward_len = int(lines[9]) - int(lines[8])
        reverse_len = int(lines[11]) - int(lines[10])

        if amplicon_chrom not in amplicon_dict: 
            amplicon_dict[amplicon_chrom] =[int(lines[9]), int(lines[10])]
        else:
            amplicon_dict[amplicon_chrom].append(int(lines[9]))
            amplicon_dict[amplicon_chrom].append(int(lines[10]))         

        for i in xrange(intargat_reads_depth):
            read1 = CreateRead()
            read2 = CreateRead()
            seq1 = genome_dict[amplicon_chrom][amplicon_start-1:(amplicon_start+151-1)].upper()
            primer1 = seq1[:forward_len]
            M_seq1 = noise(noise_rate, primer1)
            status1, seq1 = M_seq1[0], (M_seq1[1] + seq1[forward_len:(150-M_seq1[0])])
            if (len(M_seq1[1]) + insert_len) >= 150:
                read_start, read_end = insert_start, insert_start + 150 - len(M_seq1[1])
            else:
                read_start, read_end = insert_start, insert_end
            quality1 = "?"*150
            qname1 = "@ST:" + str(pos) + ":" + str(i) + " " + "1:N:0:NTCGAATC" + \
                     "_" + str(status1) + "_" + amplicon_chrom + "_" + str(read_start) + "_" + str(read_end)
            read1.qnameassign(qname1),read1.seqassign(seq1),read1.qualityassign(quality1)

            seq2 = rev_complement(genome_dict[amplicon_chrom][(amplicon_end-151-1):amplicon_end-1]).upper()
            primer2 = seq2[:reverse_len]
            M_seq2 = noise(noise_rate, primer2)
            status2, seq2 = M_seq2[0], M_seq2[1] + seq2[reverse_len:(150-M_seq2[0])]
            if (len(M_seq2[1])+insert_len) >= 150:
                read_start, read_end = insert_end-(150-len(M_seq2[1])), insert_end
            else:
                read_start, read_end = insert_start, insert_end
            qname2 = "@ST:" + str(pos) + ":" + str(i) + " " + "2:N:0:NTCGAATC" + \
                     "_" + str(status2) + "_" + amplicon_chrom + "_" + str(read_start) + "_" + str(read_end)
            quality2 = "?"*150
            read2.qnameassign(qname2),read2.seqassign(seq2),read2.qualityassign(quality2)
            Fp.write(read1, read2)

    for key,value in amplicon_dict.items():
        pos += 1
        for n in xrange(1,6,1):
            for j in xrange(oftarget_reads_depth):
                value_max = max(value)
                read1 = CreateRead()
                read2 = CreateRead()
                qname1 = "@ST-Offtarget:" + str(pos) + ":"  + str(j) + " " + "1:N:0:NTCGAATC" + "_" + "O"
                seq1 = (genome_dict[key][(value_max+500*n):(value_max+500*n+150)]).upper()
                quality1 = "?"*150
                read1.qnameassign(qname1),read1.seqassign(seq1),read1.qualityassign(quality1)
                qname2 = "@ST-Offtarget:" + str(pos) + ":"  + str(j) + " " + "2:N:0:NTCGAATC" + "_" + "O"
                seq2 = (rev_complement(genome_dict[key][(value_max+500*n+150):(value_max+500*n+300)])).upper()
                quality2 = "?"*150
                read2.qnameassign(qname2),read2.seqassign(seq2),read2.qualityassign(quality2)
                Fp.write(read1, read2)
        
if __name__ == '__main__':
    main()
