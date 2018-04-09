#!/usr/bin/python

import sys

if len(sys.argv) < 3:
	print("fastg2fasta.py inputfile outfile [minlen]")
	print("note: default minlen is set to 0; need to set it to 200, e.g., to reproduce megahit final.contigs.fa")
	sys.exit()

input_file = sys.argv[1]
output_file = sys.argv[2]

minlen = 0 #default minimum length of contigs to output; set it to 200 to reproduce metahit final.contigs.fa! 
if len(sys.argv) > 3:
	minlen = int(sys.argv[3])

ow = open(output_file,'w')
inf = open(input_file, 'r')
valid = False
seqid, seqseq = "", ""
for aline in inf:
	if aline[0] == '>':
		if valid and len(seqseq) >= minlen:
			ow.write('>' + seqid + '\n' + seqseq + '\n')
		valid = not valid
		seqid = aline[1:-1].split(':')[0].replace(';','')
		seqseq = ""
	else:
		seqseq += aline.strip()
if valid and length(seqseq) >= minlen:
	ow.write('>' + seqid + '\n' + seqseq + '\n')
ow.close()
inf.close()
