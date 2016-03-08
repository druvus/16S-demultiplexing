#!/usr/bin/env python

from Bio import SeqIO
from collections import defaultdict
import os, sys

fastqfile=sys.argv[1]
indexfile=sys.argv[2]
basename=os.path.splitext(fastqfile)[0]

d = defaultdict(dict)
i = defaultdict(dict)

print 'Loading fastq file: ' + fastqfile
handle = open(fastqfile, "rU")
for record in SeqIO.parse(handle, "fastq") :
    barcode = record.description.split(' ')[3].split('=')[1]
    #print record.description.split(' ')

    d[barcode][record.id] = record
    

handle.close()

print 'Loading index file: ' + indexfile
handle = open(indexfile, "rU")
for record in SeqIO.parse(handle, "fastq") :
    #print record.id
    i[record.id] = record

    

handle.close()

print '\nBarcodes in fastq file:\n%s' % (" ".join(d.keys()))

for key in d:
    outfile=basename + "_" + key + ".fastq"
    output_handle = open(outfile, "w")
    outIndex=basename + "_" + key + "_barcode.fastq"
    outIndex_handle = open(outIndex, "w")
    for sequence in d[key]:
        seqid=d[key][sequence].description.split(' ')[1]
        #print seqid
        
        d[key][sequence].id=seqid
        d[key][sequence].description=''
        i[seqid].description=''
        SeqIO.write(i[seqid], outIndex_handle, "fastq")
        SeqIO.write(d[key][sequence], output_handle, "fastq")
    output_handle.close()
    outIndex_handle.close()
    
#output_handle = open("example.fastq", "w")
#SeqIO.write(sequences, output_handle, "fastq")
#output_handle.close()
