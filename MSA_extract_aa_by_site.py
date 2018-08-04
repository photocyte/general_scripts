import argparse
import Bio
import Bio.SeqIO
import re

parser = argparse.ArgumentParser()
parser.add_argument('--msa',required=True,type=str,help='The coding nucleotide multiple sequence alignment')
parser.add_argument('--sites',required=True,type=str,help='A newline delimited textfile with the sites')
parser.add_argument('--seqtype',default=None,type=str,help='Either aa or nt, for nucleotide or amino acid respectively')

args = parser.parse_args()

notDNA = re.compile("[^ATCGatcgNn-]")

records = list(Bio.SeqIO.parse(args.msa,"fasta"))
sites = []
handle = open(args.sites,"rU")
for line in handle.readlines():
    sites.append(int(line.strip()))
sites.sort() ##Ensure sites are in ascending order
handle.close()

for r in records:
    if args.seqtype == None:
        results = notDNA.search(str(r.seq))
        if results != None:
            args.seqtype = "aa"
        else:
            args.seqtype = "nt"
        
    print(">"+r.name,args.seqtype)
    theSeq = Bio.Seq.Seq('')
    if args.seqtype == "nt":
        for s in sites:
            nts = r.seq[s*3-3:s*3]
            if nts == "---":
                aa = "-"
            else:
                aa=nts.translate() 
            theSeq += aa
        print(theSeq)
    elif args.seqtype == "aa":
        for s in sites:
            theSeq += r.seq[s-1:s]
        print(theSeq)    



