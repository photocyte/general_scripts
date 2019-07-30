#! /usr/bin/python 
#Timothy R. Fallon 2015

import csv
import Bio
import Bio.SearchIO
import argparse
import matplotlib.colors
import matplotlib.cm

def strip_suffix(text, suffix):
    if not text.endswith(suffix):
        return text
    else:
        return text[:len(text)-len(suffix)]

##Point of this script is to take a (complicated) blastxml file and boil it down to a simple GFF file that can be viewed in a genome browser

parser = argparse.ArgumentParser()
parser.add_argument('blast_output_xml',type=str,help='The blast output xml (input file)')

parser.add_argument('-v',
    action='store_true',
    help='Enable verbosity')

args = parser.parse_args()

if args.blast_output_xml.endswith(".gz"): 
    result_handle = gzip.open(args.blast_output_xml)
else:
    ##Try and fallback 
    result_handle = open(args.blast_output_xml)

query_results = Bio.SearchIO.parse(result_handle,'blast-xml')

output_filename = args.blast_output_xml+".gff"
handle = open(output_filename,"w")
writer = csv.writer(handle, delimiter='\t', dialect='excel')
i = 0
interval = 100
for query in query_results:
    print(query.seq_len)
    norm = matplotlib.colors.Normalize(vmin=1, vmax=query.seq_len, clip=True)
    mapper = matplotlib.cm.ScalarMappable(norm=norm, cmap=matplotlib.cm.gist_rainbow)
    for hit in query.hits:
        GFFrow = ["HITID",".","TYPE","STARTCOORD","ENDCOORD",".","STRAND",".","ID=X"]
        for hsp in hit.hsps:
            description="query_start:"+str(hsp.query_start+1)+" "+"query_end:"+str(hsp.query_end+1)
            start = 1e9
            end = 0
            color = mapper.to_rgba((hsp.query_start+1+hsp.query_end+1)/2.0,alpha=None)
            color_hex = '#'+("%0.2X" % int(color[0]*255)+("%0.2X" % int(color[1]*255))+("%0.2X" % int(color[2]*255)))
            if hsp.hit_start < start:
                start = hsp.hit_start+1
            if hsp.hit_end > end:
                end = hsp.hit_end+1
            strand = hsp.hit_strand
            GFFrow[0] = hit.id
            GFFrow[2] = "blast_hit"
            GFFrow[3] = start
            GFFrow[4] = end
            if strand == 1:
                GFFrow[6] = "+"
            if strand == -1:
                GFFrow[6] = "-"
            GFFrow[8] = "ID="+query.id+str(i)+";Description="+description+";Color="+str(color_hex)
            writer.writerow(GFFrow)
            i+=1
            if i % interval == 0:
                print(i)
                interval = interval * 10
print("Wrote to file",output_filename)
