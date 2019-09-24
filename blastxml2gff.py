#! /usr/bin/python 
#Timothy R. Fallon 2015

import csv
import Bio
import Bio.SearchIO
import argparse
import matplotlib.colors
import matplotlib.cm
import re
import random

def strip_suffix(text, suffix):
    if not text.endswith(suffix):
        return text
    else:
    	return text[:len(text)-len(suffix)]


##

def hsp_fragment_to_gff_line(y):
    q_start = y["q_start"]
    q_end = y["q_end"]
    h_start = y["h_start"]
    h_end = y["h_end"]
    strand = y["strand"]
    description="query_start:"+str(q_start)+" "+"query_end:"+str(q_end)
    color = mapper.to_rgba((q_start+q_end)/2.0,alpha=None)
    color_hex = '#'+("%0.2X" % int(color[0]*255)+("%0.2X" % int(color[1]*255))+("%0.2X" % int(color[2]*255)))
    GFFrow[0] = hit.id
    GFFrow[2] = "blast_hit"
    GFFrow[3] = h_start
    GFFrow[4] = h_end
    if strand == 1:
        GFFrow[6] = "+"
    if strand == -1:
        GFFrow[6] = "-"
    GFFrow[8] = "ID="+query.id+str(id(y))+str(random.randint(1e6,1e7))+";Description="+description+";Color="+str(color_hex)
    return GFFrow

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
handle = open(output_filename,"wb")
writer = csv.writer(handle, delimiter='\t', dialect='excel')
i = 0
interval = 100
for query in query_results:
	print query.seq_len
	norm = matplotlib.colors.Normalize(vmin=1, vmax=query.seq_len, clip=True)
        mapper = matplotlib.cm.ScalarMappable(norm=norm, cmap=matplotlib.cm.gist_rainbow)
	for hit in query.hits:
		GFFrow = ["HITID",".","TYPE","STARTCOORD","ENDCOORD",".","STRAND",".","ID=X"]
		for hsp in hit.hsps:
                        q_aln = str(hsp.aln_all[0][0].seq)
                        s_aln = str(hsp.aln_all[0][1].seq)
                        #print("_______________________________________________")
                        #print("Q:"+q_aln)
                        #print("S:"+s_aln)
                        intron_gaps = list(re.finditer("-{5,1000}",q_aln))
                        #print(intron_gaps)
                        fragments_to_make = len(intron_gaps)+1
                        #print(hsp.query_start+1)
                        #print(hsp.query_end+1)
                        #print(hsp.query_strand)
                        #print(hsp.hit_start+1)
                        #print(hsp.hit_end+1)
                        #print(hsp.hit_strand) ##If this is -1, it is a reverse oriented hit. 1 is forward oriented hit
                        #print("")
                        
                        ###Make a lookup table for the query coordinates to the nuclotide reference coordinates
                        aln_index_to_nucl_index = dict()
                        if hsp.hit_strand == 1:
                            nucl_index = hsp.hit_start
                            for aln_i in range(0,len(q_aln)):
                                aln_index_to_nucl_index[aln_i] = nucl_index
                                if s_aln[aln_i] != "-":
                                    nucl_index += 3
                                else:
                                    aln_index_to_nucl_index[aln_i] = "X"
                        elif hsp.hit_strand == -1:
                            nucl_index = hsp.hit_end
                            for aln_i in range(0,len(q_aln)):
                                aln_index_to_nucl_index[aln_i] = nucl_index
                                if s_aln[aln_i] != "-":
                                    nucl_index -= 3
                                else:
                                    aln_index_to_nucl_index[aln_i] = "X"
                        #print(aln_index_to_nucl_index)
                        #exit()
                        fragments = []
                        if fragments_to_make > 1:
                            for f in range(0,fragments_to_make):
                                if f == 0:
                                    q_index_start = 0
                                    q_index_end = intron_gaps[f].span()[0]
                                elif f > 0 and not f >= fragments_to_make-1:
                                    q_index_start = intron_gaps[f-1].span()[1]
                                    q_index_end = intron_gaps[f].span()[0]
                                elif f == fragments_to_make-1:
                                    q_index_start = intron_gaps[f-1].span()[1]
                                    q_index_end = len(q_aln)
                                #print(q_aln[q_index_start:q_index_end],q_index_start,q_index_end)
                                ##Now, have to convert the query numbering, to numbering on the nucleotide reference.
                                nucl_start = aln_index_to_nucl_index[q_index_start]
                                nucl_end = aln_index_to_nucl_index[q_index_end-1]
                                #print(s_aln[q_index_start:q_index_end],nucl_start,nucl_end,hsp.hit_start,hsp.hit_end)
                                frag = dict ()
                                frag["q_start"] = hsp.query_start+1
                                frag["q_end"] = hsp.query_end+1
                                frag["strand"] = hsp.hit_strand
                                frag["h_id"] = hit.id
                                if hsp.hit_strand == 1:
                                    frag["h_start"] = nucl_start+1
                                    frag["h_end"] = nucl_end+3
                                elif hsp.hit_strand == -1:
                                    frag["h_start"] = nucl_end-2
                                    frag["h_end"] = nucl_start-1
                                fragments.append(frag)
                        else:
                            frag = dict()
                            frag["q_start"] = hsp.query_start+1
                            frag["q_end"] = hsp.query_end+1
                            if hsp.hit_strand == 1:
                                frag["h_start"] = hsp.hit_start+1
                                frag["h_end"] = hsp.hit_end+1
                            elif hsp.hit_strand == -1:
                                frag["h_start"] = hsp.hit_start+1
                                frag["h_end"] = hsp.hit_end+1
                            frag["strand"] = hsp.hit_strand
                            frag["h_id"] = hit.id
                            fragments.append(frag)

                        for y in fragments:
                            writer.writerow(hsp_fragment_to_gff_line(y))
print "Wrote to file",output_filename
