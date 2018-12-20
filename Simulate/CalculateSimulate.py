#!/usr/bin/python2.7

import sys
import pdb
import gzip
from Bio import SeqIO
import os.path

# global variables
Result = { 'Sensitivity': 0,    # sensivity of the program
           'Specificity': 0,    # specificity of the program
           'Falsepositive': 0,  # false positive of the program
           'Falsenegative': 0,  # false negative of the program
           'Accuracy': 0,       # the accuracy of the program
           'Normal_reads': 0,   # reads with primer
           'Offtarget_reads':0  # offtarget reads
}

class FqIterator(object):
    def __init__(self, inFile):
        if inFile[-2:] == "gz":
           self.source = gzip.open(inFile, "rb")
        else:
           self.source = open(inFile, "rb")
        self.eof = False

    def next(self):
        group = []
        for i in xrange(4):
            try:
                line = self.source.next()
            except StopIteration:
                self.eof = True
                return
            group.append( line.rstrip() )

        return group

    def close(self):
        self.source.close()

def RevComplement(seq):
        seq = seq.upper()
        basecomplement = {'A':'T',
                          'C':'G',
                          'G':'C',
                          'T':'A',
                          '-':'-',
                          'N':'N'}
        letters = list(seq)
        letters = [basecomplement[base] for base in letters]
        complement = (''.join(letters))                                                 #gives the complement of the bases in list letters

        return complement[::-1] 

def load_reference(fa_file):
    fa_dict = {}

    fasta_seq = SeqIO.parse(fa_file, "fasta")

    for fasta in fasta_seq:
        key = fasta.id
        value = fasta.seq
        if key not in fa_dict:
            fa_dict[key] = value
        else:
            sys.stderr.wirte("Err: Don't have the same chromsome!!!\n")
            sys.exit(-1)

    return fa_dict


def get_seq(chrom, start, end, fa_dict):
    return fa_dict[chrom][start-1:end]


def calculate(read, fa_dict):
    ''' read = [seqname, sequence, +, quality]
    '''
    if read[0][-1] != 'O':
        Result['Normal_reads'] += 1      #normal reads with primer
        if len(read[1]) < 150: 
            Result['Sensitivity'] += 1   # has primer and detected by the program
        else:
            Result['Falsenegative'] += 1 # has primer but failed to detect

        name = read[0].rstrip().split('_')

        flag = read[0].rstrip().split(" ")[1].split(":")[0]
        if flag == "1":
            insert_seq = get_seq(name[2], int(name[3]), int(name[4]), fa_dict)
        else:
            insert_seq = RevComplement(get_seq(name[2], int(name[3]), int(name[4]), fa_dict))
        if insert_seq.upper() == read[1]:
            Result['Accuracy'] += 1
    else:
        Result['Offtarget_reads'] += 1   #offtarget reads without primer
        if len(read[1]) < 150:
            Result['Falsepositive'] += 1 # without primer but program think it has
        else:
            Result['Specificity'] += 1   # without primer and program realize that


def main():
    args = sys.argv

    if len(args) != 4:
        sys.stderr.write("Usage: python %s <ref.fa> <read1> <read2>\n" %args[0])
        sys.exit(-1)

    #pdb.set_trace()
    Fastq1 = FqIterator(args[2]) # read1 of the paired fastq file
    Fastq2 = FqIterator(args[3]) # read2 of the paired fastq file
    fa_dict = load_reference(args[1])

    while (True):
        read1 = Fastq1.next()
        read2 = Fastq2.next()
        if (not Fastq1.eof) and (not Fastq2.eof):
            calculate(read1, fa_dict); calculate(read2, fa_dict)
        else:
            break

    Fastq1.close(); Fastq2.close()


    # Summary
    sys.stdout.write("Sensitivity: %.2f\n" %(Result['Sensitivity']/float(Result['Normal_reads'])*100))
    sys.stdout.write("Specificity: %.2f\n" %(Result['Specificity']/float(Result['Offtarget_reads'])*100))
    sys.stdout.write("Falsepositive: %.2f\n" %(Result['Falsepositive']/float(Result['Offtarget_reads'])*100))
    sys.stdout.write("Falsenegative: %.2f\n" %(Result['Falsenegative']/float(Result['Normal_reads'])*100))
    sys.stdout.write("Accuracy: %.2f\n" %(Result['Accuracy']/float(Result['Normal_reads'])*100))



if __name__ == '__main__':
    main()
