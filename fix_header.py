#!/usr/bin/env python

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict
import os, sys

#################################################
def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    complseq = [complement[base] for base in seq]
    return complseq

def reverse_complement(seq):
    seq = list(seq)
    seq.reverse()
    return ''.join(complement(seq))

#################################################

mappingfile=sys.argv[1]
fastafile=sys.argv[2]
barcode=sys.argv[3]
directory_name=sys.argv[4]

basename = os.path.basename(fastafile)
outname = basename.replace('seqs', 'corrected')

d = defaultdict(dict)

print 'Loading mapping file: ' + mappingfile

handle = open(mappingfile, "rU")
header = handle.readline()
for record in handle :
    tmp = record.split('\t')
    ID = tmp[0]
    fprimer = tmp[13]
    rprimer = tmp[14]
    d[rprimer][fprimer] = ID

handle.close()


tmp = fastafile.split('/')[1]
tmp=os.path.splitext(tmp)[0]
barcode_rc = reverse_complement(barcode)


print 'Found reverse barcode: ' + barcode
print 'Reverse complemented reverse barcode: ' + barcode_rc

print 'Loading fasta file: ' + fastafile
handle = open(fastafile, "rU")
out = '%s/%s' % (directory_name, outname)
outfile = open(out, 'w')
skipline = False

for line in handle:
    if line.startswith( '>' ):
        ID, machine, orig_bc, new_bc, bc_diffs = line.split(' ')
        fprimer = reverse_complement(new_bc.split('=')[1])
        rprimer = barcode_rc
        try:
            ID = d[rprimer][fprimer]
        except KeyError:
            skipline = True
            continue
        bc = rprimer + fprimer
        headerline = '>' + ID + ' ' + machine + ' orig_bc=' + bc + ' new_bc=' + bc + ' bc_diffs=0\n'
        outfile.write(headerline)
    else:
        if skipline != True:
            outfile.write(line)
        skipline = False
      
      

handle.close()
outfile.close()

print 'Sequences are written to: ' + outname + '\n'
