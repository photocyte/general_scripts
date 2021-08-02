from Bio import SeqIO
from Bio import GenBank
import os
import argparse

## Keywords
# biopython genbank to peptide fasta 
# Updated from: https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank2fasta/

parser = argparse.ArgumentParser()
parser.add_argument("genbank_filename")
parser.add_argument("--type",default="CDS")
parser.add_argument("--no_filter",default=False, action='store_true')
args = parser.parse_args()

gbk_filename = args.genbank_filename
gff_filename = gbk_filename+"_converted.gff"

input_handle  = open(gbk_filename, "r")
output_handle = open(gff_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print("Dealing with GenBank record %s" % seq_record.id)
    for seq_feature in seq_record.features :
        if seq_feature.type==args.type or args.no_filter:
            print(seq_feature)
            qualifiers = [ str(k)+"="+str(v) for k,v in seq_feature.qualifiers.items() ]
            qualifiers.sort()
            qualifiers_str = ";".join(qualifiers)
            if seq_feature.strand == 1:
                strand = "+"
            elif seq_feature.strand == -1:
                strand = "+"
            else:
                strand = "?"
            gff_fields = [seq_record.id,"unknown",args.type,seq_feature.location.start,seq_feature.location.end,".",strand,".",qualifiers_str]
            gff_fields = [str(x) for x in gff_fields]
            output_handle.write("\t".join(gff_fields)+os.linesep)
            #assert len(seq_feature.qualifiers['translation'])==1
        else:
            pass
            #print(seq_feature)            
#            output_handle.write(">%s from %s\n%s\n" % (
#                   seq_feature.qualifiers['locus_tag'][0],
#                   seq_record.name,
#                   seq_feature.qualifiers['translation'][0]))

output_handle.close()
input_handle.close()
print("Done.")
