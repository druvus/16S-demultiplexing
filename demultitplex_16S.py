#!/usr/bin/env python

import subprocess
import sys
import os
import argparse

parser = argparse.ArgumentParser()
parser._optionals.title = "Options"
parser.add_argument('-m', nargs = 1, help = 'Mapping file. [Required]', required=True)
parser.add_argument('-i', nargs = 1, help = 'Path to the R1, R2 and I1 and I2 files. [Required]', required=True)
parser.add_argument('-o', nargs = 1, help = 'A name that is given to the output files. [Required]', required=True)
parser.add_argument('--no_otu_picking', action ='store_true', help = 'Add if no OTU-picking should be performed.')
parser.add_argument('--keep_temp', action ='store_true', help = 'If given, the files that are produced during the script will be kept.')
args = parser.parse_args()
arguments = vars(args)

script_dir = "/scratch/andrsjod/git_code/16S-demultiplexing"
shortname = "ratdata"

r1 = '%s/%s_R1.fastq.tar.gz' % (''.join(arguments['i']), shortname)
r2 = '%s/%s_R2.fastq.tar.gz' % (''.join(arguments['i']), shortname)
i1 = '%s/%s_I1.fastq.tar.gz' % (''.join(arguments['i']), shortname)
i2 = '%s/%s_I2.fastq.tar.gz' % (''.join(arguments['i']), shortname)
sample = ''.join(arguments['o'])
directory_name = ''.join(arguments['i']).split('/')[-2]


if arguments['no_otu_picking']:
    print(""" ********************************************************

              Closed OTU-picking will not be performed in the script.

              ********************************************************""")
    otu = False
else:
    otu=True
if arguments['keep_temp']:
    print(""" ********************************************************

              All files created in the script will be kept in the foler deplex_temporary_files.

              ********************************************************""")
    keep=True
else:
    keep = False

I1 = os.path.basename(i1).rsplit('.', 2)[0]
I2 = os.path.basename(i2).rsplit('.', 2)[0]
R1 = os.path.basename(r1).rsplit('.', 2)[0]
R2 = os.path.basename(r2).rsplit('.', 2)[0]

def shell_source(script):
    pipe = subprocess.Popen(". %s; env" % script, stdout=subprocess.PIPE, shell=True)
    output = pipe.communicate()[0]
    env = dict((line.split("=", 1) for line in output.splitlines()))
    os.environ.update(env)

print("Demultiplexing started...")
#source = "/home/andrsjod/qiime_software_1.7/activate.sh"
#shell_source(source);

cmd = "mkdir deplex_temporary_files"
subprocess.call(cmd, shell=True)

## Unpack and move index files to current folder
cmd1 = "tar  --to-stdout -xvf %s > deplex_temporary_files/%s" % (i1, I1)
subprocess.call(cmd1, shell=True)
cmd2 = "tar  --to-stdout -xvf %s > deplex_temporary_files/%s" % (i2, I2)
subprocess.call(cmd2, shell=True)

## Unpack and move read files to current folder
cmd1 = "tar  --to-stdout -xvf %s > deplex_temporary_files/%s" % (r1, R1)
subprocess.call(cmd1, shell=True)
cmd2 = "tar  --to-stdout -xvf %s > deplex_temporary_files/%s" % (r2, R2)
subprocess.call(cmd2, shell=True)


cmd = "validate_mapping_file.py -m %s -o deplex_temporary_files/check_map" % ''.join(arguments['m'])
subprocess.call(cmd, shell=True)

cmd = "%s/fix_mappingfile.py deplex_temporary_files/check_map/*_corrected.txt deplex_temporary_files" % script_dir
subprocess.call(cmd, shell=True)

## Demultiplexing and get quality filtering.
cmd1 = "split_libraries_fastq.py -o deplex_temporary_files/mapping_%s_1 -i deplex_temporary_files/%s -b deplex_temporary_files/%s   --rev_comp_mapping_barcodes -m deplex_temporary_files/*_corrected_1.txt  --store_demultiplexed_fastq" % (sample, R1, I1)
cmd2 = "split_libraries_fastq.py -o deplex_temporary_files/mapping_%s_2 -i deplex_temporary_files/%s -b deplex_temporary_files/%s   --rev_comp_mapping_barcodes -m deplex_temporary_files/*_corrected_1.txt  --store_demultiplexed_fastq" % (sample, R2, I1)

print("In progress: demultiplexing and quality filtering of the R1 data...")
subprocess.call(cmd1, shell=True)
print("In progress: demultiplexing and quality filtering of the R2 data...")
subprocess.call(cmd2, shell=True)


## Syncing
cmd = "%s/syncsort_fq deplex_temporary_files/mapping_%s_1/seqs.fastq deplex_temporary_files/mapping_%s_2/seqs.fastq deplex_temporary_files" % (script_dir, sample, sample)
subprocess.call(cmd, shell=True)

## Cutadapt
## 16S ATTAGAWACCCBDGTAGTCC TTACCGCGGCKGCTGGCAC Since there are variation among the primers it can be helpful to set the -e (allowed error rate) parameter a little less stringent. default = 0.1
cmd = "cutadapt -a ATTAGATACCCTAGTAGTCC deplex_temporary_files/mapping_%s_1/seqs.fastq_paired.fq -o deplex_temporary_files/%s_only_R1_clean.fq -e 0.2" % (sample, sample)
subprocess.call(cmd, shell=True)
cmd = "cutadapt -a TTACCGCGGCTGCTGTCAC deplex_temporary_files/mapping_%s_2/seqs.fastq_paired.fq -o deplex_temporary_files/%s_only_R2_clean.fq -e 0.2" % (sample, sample)
subprocess.call(cmd, shell=True)

## Flash
cmd = "flash -m 10 -M 250 -r 250 -f 253 -t 12 -o deplex_temporary_files/%s deplex_temporary_files/%s_only_R1_clean.fq  deplex_temporary_files/%s_only_R2_clean.fq" % (sample, sample, sample)
subprocess.call(cmd, shell=True)

## Syncing
cmd = "%s/syncsort_fq_readindex deplex_temporary_files/%s.extendedFrags.fastq deplex_temporary_files/%s deplex_temporary_files/%s" % (script_dir, sample, I1, I2)
subprocess.call(cmd, shell=True)

## Splitting
print("In progress: split_libraries_fastq.py...")
cmd = "split_libraries_fastq.py -o deplex_temporary_files/mapping_%s -i deplex_temporary_files/%s.extendedFrags.fastq_synced.fq -b deplex_temporary_files/%s_synced.fq  --rev_comp_mapping_barcodes -m deplex_temporary_files/*_corrected_1.txt  --store_demultiplexed_fastq --phred_offset 33 -p 0.6" % (sample, sample, I1)
subprocess.call(cmd, shell=True)

cmd = "%s/split_fastq.py deplex_temporary_files/mapping_%s/seqs.fastq deplex_temporary_files/%s_synced.fq" % (script_dir, sample, I2)
out = subprocess.Popen(cmd, stdout = subprocess.PIPE, shell=True).communicate()
barcodes = out[0].split('\n')[4].split(' ')

for barcode in barcodes:
    print ('Analysing ' + barcode)
    cmd = "split_libraries_fastq.py -o deplex_temporary_files/mapping_%s -i deplex_temporary_files/mapping_%s/seqs_%s.fastq -b deplex_temporary_files/mapping_%s/seqs_%s_barcode.fastq   --rev_comp_mapping_barcodes -m deplex_temporary_files/*_corrected_2.txt  --store_demultiplexed_fastq --rev_comp_barcode --phred_offset 33 -p 0.6" % (barcode, sample, barcode, sample, barcode)
    subprocess.call(cmd, shell=True)
    cmd = "chmod 755 -R deplex_temporary_files/mapping_%s" % barcode
    subprocess.call(cmd, shell=True)

cmd = "mkdir deplex_temporary_files/mapping_all"
subprocess.call(cmd, shell=True)
cmd = "mkdir deplex_temporary_files/mapping_all_fastq"
subprocess.call(cmd, shell=True)

for barcode in barcodes:
    print ('Copying ' + barcode)
    cmd = "cp deplex_temporary_files/mapping_%s/seqs.fna deplex_temporary_files/mapping_all/seqs_%s.fna" % (barcode, barcode)
    subprocess.call(cmd, shell=True)
    cmd = "cp deplex_temporary_files/mapping_%s/seqs.fastq deplex_temporary_files/mapping_all_fastq/seqs_%s.fastq" % (barcode, barcode)
    subprocess.call(cmd, shell=True)

for barcode in barcodes:
    print ('Fixing fasta and fastq file for ' + barcode)
    cmd = "%s/fix_header.py deplex_temporary_files/check_map/*_corrected.txt deplex_temporary_files/mapping_all/seqs_%s.fna %s deplex_temporary_files" % (script_dir, barcode, barcode)
    subprocess.call(cmd, shell=True)
    cmd = "%s/fix_header_fastq.py deplex_temporary_files/check_map/*_corrected.txt deplex_temporary_files/mapping_all_fastq/seqs_%s.fastq %s deplex_temporary_files" % (script_dir, barcode, barcode)
    subprocess.call(cmd, shell=True)



## Fix final fasta and fastq files
cmd = "cat deplex_temporary_files/corrected_*.fna > deplex_temporary_files/corrected_tmp.fna"
subprocess.call(cmd, shell=True)
cmd = "cat deplex_temporary_files/corrected_*.fastq > deplex_temporary_files/corrected_tmp.fastq"
subprocess.call(cmd, shell=True)
cmd = "%s/fix_ID.py deplex_temporary_files/corrected_tmp.fna deplex_temporary_files" % (script_dir)
subprocess.call(cmd, shell=True)
cmd = "%s/fix_ID_fastq.py deplex_temporary_files/corrected_tmp.fastq deplex_temporary_files" % (script_dir)
subprocess.call(cmd, shell=True)

## TODO

cmd = "mkdir processed_data"
subprocess.call(cmd, shell=True)

## Move the created files that should be kept to new locations.
cmd = "mkdir processed_data/%s" % directory_name
subprocess.call(cmd, shell=True)
cmd = "mkdir processed_data/%s/16S" % directory_name
subprocess.call(cmd, shell=True)
cmd = "mv deplex_temporary_files/corrected_all.* processed_data/%s/16S/." % directory_name
subprocess.call(cmd, shell=True)

cmd = "mkdir %s_demultiplexed_data" % sample
subprocess.call(cmd, shell=True)
cmd = "mkdir %s_demultiplexed_data/Histograms" % sample
subprocess.call(cmd, shell=True)
cmd = "mv deplex_temporary_files/%s.hist* %s_demultiplexed_data/Histograms/." % (sample, sample) 
subprocess.call(cmd, shell=True)
cmd = "mv deplex_temporary_files/check_map/*_corrected.txt %s_demultiplexed_data/." % sample
subprocess.call(cmd, shell=True)
cmd = "cp %s_demultiplexed_data/*_corrected.txt processed_data/%s/16S/." % (sample, directory_name)
subprocess.call(cmd, shell=True)

## Remove the files that are created during the script. 
if keep == False:
    cmd = "rm -R deplex_temporary_files"
    subprocess.call(cmd, shell=True)

## Create links to the fasta file and fastq file.
cmd = "ln -sf  ../processed_data/%s/16S/corrected_all.fna %s_demultiplexed_data/%s.fna" % (directory_name, sample, sample)
subprocess.call(cmd, shell=True)
cmd = "ln -sf  ../processed_data/%s/16S/corrected_all.fastq %s_demultiplexed_data/%s.fastq" % (directory_name, sample, sample)
subprocess.call(cmd, shell=True)

cmd = "cat %s_demultiplexed_data/*_corrected.txt | cut -f 1" % sample
complete = subprocess.Popen(cmd, stdout = subprocess.PIPE, shell=True).communicate()
all_samples = filter(None, complete[0].split('\n')[1:])

f = open('%s_demultiplexed_data/%s_number_of_reads.txt' % (sample, sample), 'w')
f.write('%-20s\t\t%10s\n' % ('#SampleID', 'Number of reads'))
for sampleID in all_samples:
    cmd = "grep %s %s_demultiplexed_data/%s.fna | wc -l" % (sampleID, sample, sample)
    line = subprocess.Popen(cmd, stdout = subprocess.PIPE, shell=True).communicate()
    f.write('%-20s \t\t %10s\n' % (sampleID, line[0].split('\n')[0]))

f.close()

print(""" ********************************************************

           Demultiplexing is done, the fasta file is now available!

           ********************************************************""")

## Closed OTU-picking
if otu == True:
    print ("Closed OTU-picking is initiated...")
    cmd = "pick_closed_reference_otus.py -i %s_demultiplexed_data/%s.fna -o %s_demultiplexed_data/Closed_OTU_picking -r /scratch/andrsjod/qiime_db/gg_13_8_otus/rep_set/97_otus.fasta -t /scratch/andrsjod/qiime_db/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt -aO 12" % (sample, sample, sample)
    subprocess.call(cmd, shell=True)
'''
    ## Filter singletons
    print "Closed OTU-picking is completed, filtration of singletons is initiated..."
    cmd = "filter_otus_from_otu_table.py -i %s_demultiplexed_data/Closed_OTU_picking/otu_table.biom -o %s_demultiplexed_data/%s_n2.biom -n 2" % (sample, sample, sample)
    subprocess.call(cmd, shell=True)

## Copy the mapping file, OTU-table and OTU-mapping file to the results folder on Sourcetracker.
cmd = "mkdir /mnt/sourcetracker/results/16S/%s" % directory_name
subprocess.call(cmd, shell=True)
cmd = "cp %s_demultiplexed_data/Closed_OTU_picking/otu_table.biom /mnt/sourcetracker/results/16S/%s/." % (sample, directory_name)
subprocess.call(cmd, shell=True)
cmd = "cp %s_demultiplexed_data/Closed_OTU_picking/uclust_ref_picked_otus/%s_otus.txt /mnt/sourcetracker/results/16S/%s/." % (sample, sample, directory_name)
subprocess.call(cmd, shell=True)
cmd = "cp %s_demultiplexed_data/*_corrected.txt /mnt/sourcetracker/results/16S/%s/." % (sample, directory_name)
subprocess.call(cmd, shell=True)
'''

print("********************************************************\n\nThe script is done, the OTU-table is now available!\n\n********************************************************")
