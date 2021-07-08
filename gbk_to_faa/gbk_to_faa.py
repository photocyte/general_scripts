from Bio import SeqIO
from Bio import GenBank
import argparse

## Keywords
# biopython genbank to peptide fasta 
# Updated from: https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/genbank2fasta/

parser = argparse.ArgumentParser()
parser.add_argument("genbank_filename")
args = parser.parse_args()

gbk_filename = args.genbank_filename
faa_filename = gbk_filename+"_converted.faa"

input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print("Dealing with GenBank record %s" % seq_record.id)
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" :
            assert len(seq_feature.qualifiers['translation'])==1
            output_handle.write(">%s from %s\n%s\n" % (
                   seq_feature.qualifiers['locus_tag'][0],
                   seq_record.name,
                   seq_feature.qualifiers['translation'][0]))

output_handle.close()
input_handle.close()
print("Done.")
