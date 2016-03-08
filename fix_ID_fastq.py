#!/usr/bin/env python
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict
import os, sys


fastqfile=sys.argv[1]
directory= sys.argv[2]
basename = os.path.basename(fastqfile)
outname = basename.replace('tmp', 'all')


print 'Loading fastq file: ' + fastqfile
handle = open(fastqfile, "rU")
out = '%s/%s' % (directory, outname)
outfile = open(out, 'w')
count = 0
skipline = False

for line in handle:
    if line.startswith( '@' ):
        try:
		ID, machine, orig_bc, new_bc, bc_diffs = line.split(' ')
		ID = ID + '_' + str(count)
        	headerline = ID + ' ' + machine + ' ' + orig_bc + ' ' + new_bc + ' ' + bc_diffs
        	outfile.write(headerline)
        	count += 1
	except:
		pass
    else:
        outfile.write(line)

handle.close()
outfile.close()

print 'Sequences are written to: ' + outname + '\n'
