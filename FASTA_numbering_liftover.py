#! /Users/tim/anaconda_ete/bin/python
import Bio
import Bio.SeqIO
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('-A', required=True,help='A FASTA record where you know the numbering')
parser.add_argument('-B', required=True,help='A FASTA record where you want to liftover the numbering to (e.g. it is gapped)')
parser.add_argument('--seqtype',default=None,help='Either aa or nt, for nucleotide or amino acid respectively')
parser.add_argument('--sites',default=None)
args = parser.parse_args()

if args.sites != None:
    sites = []
    handle = open(args.sites,"rU")
    for line in handle.readlines():
        sites.append(int(line.strip()))
    sites.sort() ##Ensure sites are in ascending order
    handle.close()

record_A = list(Bio.SeqIO.parse(args.A,"fasta"))[0]
record_B = list(Bio.SeqIO.parse(args.B,"fasta"))[0]

notDNA = re.compile("[^ATCGatcgNn-]")

if notDNA.search(str(record_A.seq)) != None or notDNA.search(str(record_B.seq)) != None:
    args.seqtype = "aa"
else:
    args.seqtype = "nt"

i=1
j=1
while i-1 < len(record_A.seq):
    while record_A.seq[i-1] == "-":
        i+=1
    while record_B.seq[j-1] == "-":
        j+=1
    if args.seqtype == "nt":
         print(record_A.seq[i-1],i/3.0,j/3.0,record_B.seq[j-1])
    if args.seqtype == "aa" and args.sites == None:
         print(record_A.seq[i-1],i,j,record_B.seq[j-1])
    elif args.sites != None and i in sites:
         print(record_A.seq[i-1],i,j,record_B.seq[j-1]) ##Only print if it is in the sites list
    assert record_A.seq[i-1] == record_B.seq[j-1]
    i+=1
    j+=1   
