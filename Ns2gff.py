#! /usr/bin/env python
#Timothy R. Fallon 2015

import csv
import Bio
import Bio.SeqIO
import argparse
import os
import os.path
import sys

##Point of this script is to take the Ns from a scaffolded fasta file, and produce a GFF file that can be viewed in a genome browser

parser = argparse.ArgumentParser()
parser.add_argument('fasta_input',type=str,help='The fasta file (input file)')

parser.add_argument('-v',
    action='store_true',
    help='Enable verbosity')

args = parser.parse_args()

if args.fasta_input.endswith(".gz"): 
    ##Decompress gzip compression on the fly
    result_handle = gzip.open(args.fasta_input)
    sys.stderr.write("Opening "+args.fasta_input+" with gzip\n")

else:
    ##Try and fallback if its not gzip compressed
    result_handle = open(args.fasta_input)
    sys.stderr.write("Opening "+args.fasta_input+'\n')

fasta_records = Bio.SeqIO.parse(result_handle,'fasta')
output_filename = "./Ns2gff_"+os.path.basename(args.fasta_input)+".gff"
handle = open(output_filename,"wb")
writer = csv.writer(handle, delimiter='\t', dialect='excel')

i = 0
interval = 10
features = []
feature = dict()

def feature_to_gffrow(feat):
    GFFrow = ["HITID",".","TYPE","STARTCOORD","ENDCOORD",".","STRAND",".","ID=X"]
    GFFrow[0] = feat['id']
    GFFrow[2] = "unknown bases"
    GFFrow[3] = str(feat['start'])
    GFFrow[4] = str(feat['stop'])
    GFFrow[6] = "+"
    GFFrow[8] = "ID="+feat['id']+"-"+feat['count']
    return GFFrow

sys.stderr.write("Collecting Ns from FASTA records and writing to standard out...\n")
count = 0
for record in fasta_records:
    feature['status'] = False
    feature['start'] = '' 
    feature['stop'] = ''
    feature['id'] = ''
    feature_count = 0
    ##Start scanning record.seq for Ns
    for i in range(0,len(record.seq)):
        if feature['status'] == False and (record.seq[i] == "N" or record.seq[i] == "n"):
            feature['status'] = True
            feature['start'] = i + 1 ##GFF 1 indexed
        if feature['status'] == True and (record.seq[i] != "N" and record.seq[i] != "n"):
            feature_count += 1
            feature['status'] = 'endfound'
            feature['stop'] = i ##Should be i - 1, but GFF 1 indexed
            feature['id'] = record.id
            feature['count'] = "N"+str(feature_count)
            gffrow = feature_to_gffrow(feature)
            sys.stdout.write("\t".join(gffrow)+"\n")
            features.append(feature)
            feature = dict()
            feature['status'] = False
    count +=1
    if count % interval == 0:
        sys.stderr.write("Completed: "+str(count)+" scaffolds.\n")
        interval = interval * 10
sys.stderr.write("Found "+str(len(features))+" strings of Ns in "+str(count)+" scaffolds.\n")

exit()

