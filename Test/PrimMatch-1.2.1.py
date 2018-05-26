#!/usr/bin/python2.7

'''
    Date: 2017-11-04

    NOTE: we have modify the rules of sequence query.
    Assuming that every fastq sequence has a primer seq, we are sure
    that after query for "3+mismatch" bases, there is no possibility 
    to find the primer seq anymore!
'''

import sys
import os
import gzip
import pdb
import os.path

# global variable
MAXPRIMLEN = 50
READLEN = 150
BLEN = 6 # bufer length for search

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
            if i == 0:
                group.append( line.split()[0] )
            else:
                group.append( line.rstrip() )
        return group

    def close(self):
        self.source.close()


def GetPrim(primerfile):
    ''' primer = [('CCCTTCCACATA','CAAGTGAGGTCCTCAAAT', 99), (), ...]
    '''
    primer = []

    with open(primerfile, "r") as fp:
        while (1):
            line = fp.readline()
            if not line.startswith('#'):
                break

        for line in fp.readlines():
            llist = line.split()
            primer.append( (llist[0], llist[1], int(llist[2])) )

    return primer


def SeqIndex(primseq, loc1, loc2, kmer):
    ''' primseq: primer sequence [str]
        loc1: primseq index in the primlist [int]
        loc2: primseq index in the primer pair tuple [int]
        i: kmer seq index in the primer sequence
        index = {'ATGCG':[(loc1, loc2, i),(), ...], ...}
    '''
    index = {}
    cutlen = len(primseq) - kmer + 1

    for i in xrange(cutlen):
        key = primseq[i:i+kmer]
        if key not in index:
            index[key] = [(loc1, loc2, i)]
        else:
            index[key].append((loc1, loc2, i))

    return index


def Index(primlist, kmer=8):
    ''' index = {'ATGCG':[(0,1,0),(0,1,1)], 'TGCGA':[(),(),(),...], ...}
    '''
    index = {}
    pcount = len(primlist)

    for pindex in xrange(pcount):
        for i in [0, 1]:
            tmpindex = SeqIndex(primlist[pindex][i], pindex, i, kmer)
            for key in tmpindex:
                if key not in index:
                    index[key] = tmpindex[key]
                else:
                    for ele in tmpindex[key]:
                        index[key].append(ele)
    return index


def RevComp(sequence):
    base = {'A':'T','T':'A','G':'C','C':'G'}
    compseq = ''

    for b in sequence:
        compseq += base[b]

    return compseq[::-1]


def MisCheck(str1, str2, mismatch):
    misnum = 0
    len1 = len(str1)
    len2 = len(str2)

    minlen = len1 if len1 < len2 else len2
    for i in xrange(minlen):
        if str1[i] != str2[i]:
            misnum += 1

        if misnum > mismatch:
            return False

    return True


def SeqQuery(seq, index, primlist, mismatch, kmer):
    '''CORE FUNCTION
          seq: AGAAATTTGCGGAGTAAGTTGCGCTGGGGCTTTCGGCGGCGGCGATTTCGCC
       primer:      TTTGCGGAGTAAG
                    |           |
                  pstart       pend
       condition as above:  i = 5; h[2] = 0
    '''
    cutlen = 2 * BLEN

    for i in xrange(cutlen):
        try:
            hitlist = index[seq[i:i+kmer]]
        except KeyError:
            continue

        for h in hitlist:
            if i - h[2] >=0:
                string, primer = seq[i-h[2]:], primlist[h[0]][h[1]]
            elif (i - h[2] < 0) and (i - h[2] >= -mismatch):
                string, primer = seq, primlist[h[0]][h[1]][h[2]-i:]
            else:
                continue

            if MisCheck(string, primer, mismatch):
                pstart = 0 if i-h[2] < 0 else i-h[2]
                pend = pstart + len(primer) - 1
                return (h[0], h[1], pstart, pend)
    return ()


def PrimQuery(readseq, index, primlist, mismatch=3, kmer=8):
    fwdseq = readseq[:MAXPRIMLEN+BLEN]

    loc = SeqQuery(fwdseq, index, primlist, mismatch, kmer)
    if not loc:
        return ()

    insertlen = primlist[loc[0]][2]
    # for the condition of target sequencing
    if loc[3]+insertlen > READLEN:
        return (loc[0], loc[1], loc[3]+1, 0) # the seqend never exceed insertlen

    revstart = loc[3] + insertlen - 3
    revseq = readseq[revstart:]
    revprim = RevComp(primlist[loc[0]][loc[1]^1])

    # Indexing for the reverse primer
    revindex = SeqIndex(revprim, 0, 0, kmer)

    revloc = SeqQuery(revseq, revindex, [(revprim, 'dummy', 0)], mismatch, kmer)
    if not revloc:
        return ()

    insertstart = loc[3] + 1
    insertend = revstart + revloc[2] - 1

    return (loc[0], loc[1], insertstart, insertend)


def PrintMatch(target, Primer, seq):
    fw = Primer[target[0]][target[1]]
    rv = RevComp(Primer[target[0]][target[1]^1])

    seqshift = len(fw) - target[2]
    sys.stdout.write('%s%s\n' %(' '*seqshift, seq))
    if target[3] != 0:
        sys.stdout.write( '%s%s%s\n' %(fw,' '*(target[3]-target[2]+1), rv) )
    else:
        sys.stdout.write( '%s\n' %fw)


def main():
    args = sys.argv
    # args[1]: amplicon.txt
    # args[2]: fastq file

    if len(args) != 3:
        sys.stderr.write('Usage: python verify.py <amplicon.txt> <*.fq/*.fq.gz>\n')
        sys.exit(1)

    primer = GetPrim(args[1])
    index = Index(primer, 8)
    fastq = FqIterator(args[2])

    while ( True ):
        read = fastq.next()
        if fastq.eof:
            break

        target = PrimQuery(read[1], index, primer, 3, 8)
        if target:
            PrintMatch(target, primer, read[1])

    fastq.close


if __name__ == '__main__':
    main()
