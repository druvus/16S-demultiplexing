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
output_directory = sys.argv[2]

basename = os.path.basename(mappingfile)
outname1 = basename.replace('.txt', '_1.txt')
outname2 = basename.replace('.txt', '_2.txt')

step1 = defaultdict(dict)
step2 = defaultdict(dict)

print 'Loading mapping file: ' + mappingfile

handle = open(mappingfile, "rU")
header = handle.readline()
for record in handle :
    tmp = record.split('\t')
    step1[tmp[14]] = tmp[11] + tmp[7]
    step2[tmp[13]] = tmp[12] + tmp[8]

handle.close()
out1 = '%s/%s' % (output_directory, outname1)
out2 = '%s/%s' % (output_directory, outname2)
outfile1 = open(out1, 'w')
outfile2 = open(out2, 'w')

print 'Writing mapping file: ' + outname1

outfile1.write('#SampleID' + '\t' + 'BarcodeSequence' + '\t' + 'LinkerPrimerSequence' + '\t' + 'SampleType' + '\t' + 'Region' + '\t' + 'Ordering' + '\t' + 'Description'+ '\n')
counter = 0
for key in step1.keys():
    counter += 1
    outfile1.write('Reversebarcode' + str(counter) + '\t' + key + '\t' + step1[key] + '\t' + 'NA' + '\t' + 'NA' + '\t' + str(counter) + '\t' + 'NA'+ '\n')
    
print 'Writing mapping file: ' + outname2

outfile2.write('#SampleID' + '\t' + 'BarcodeSequence' + '\t' + 'LinkerPrimerSequence' + '\t' + 'SampleType' + '\t' + 'Region' + '\t' + 'Ordering' + '\t' + 'Description'+ '\n')
counter = 0
for key in step2.keys():
    counter += 1
    outfile2.write('Forwardbarcode' + str(counter) + '\t' + key + '\t' + step2[key] + '\t' + 'NA' + '\t' + 'NA' + '\t' + str(counter) + '\t' + 'NA'+ '\n')
    
outfile1.close()
outfile2.close()
