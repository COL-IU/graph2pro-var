#!/usr/bin/env python
#Nov 8 2018, Ye
import sys
import re

if len(sys.argv) < 6:
	print ("Usage: " + sys.argv[0] + " <fgs-ms-ide-result> <fgs-fasta-file> <graph2pep-ms-ide-result> <graph2pep-fasta-file> <output-file>")
	sys.exit()
fgsms, fgsfa, dbpepms, dbpepfa, outfile = sys.argv[1:6]
msfiles = [dbpepms, fgsms]

cidx, eidx, pidx = 0, 1, 2
des = ["DBGraph2Pro", "FragGeneScan"]
#read peptide identification results
pepHash = [{}, {}]
for i in range(2):
	f = open(msfiles[i], "r")
	for line in f:
		if line[0] == '#':
			continue
		subs = line.split('\t')
		peptide = re.sub('[^a-zA-Z]', '', subs[9])
		protein = subs[10].split('(')[0]
		if protein.startswith('XXX_'):
			continue
		evalue = float(subs[14]) 
		if peptide in pepHash[i]:
			update = []
			exists = False
			#protein group (multiple proteins contain the same peptide)
			for one in pepHash[i][peptide]:
				onenew = one[:]
				if protein == one[pidx]:
					onenew[cidx] += 1
					if evalue < onenew[eidx]:	
						onenew[eidx] = evalue
					exists = True
				update.append(onenew)
			if not exists:
				onenew = [1, evalue, protein]
				update.append(onenew)
			pepHash[i][peptide] = update
		else:
			pepHash[i][peptide] = [[1, evalue, protein],]
	f.close()

fw = open(outfile,'w')
fw.write('#peptide\tprotein\tbestEvalue\tsource\tpeptideCounts\tedge1\tstart\tedge2\tend\tstopCodon\n')

#get pep location in graph
graphinfo = {}
for pep in pepHash[0]:
	for one in pepHash[0][pep]:
		protein = one[pidx]
		graphinfo[protein] = ""
f = open(dbpepfa, "r")
name = ""
for aline in f:
	if aline[0] == '>':
		subs = aline[1:-1].split()	
		name = subs[0]
		if name in graphinfo:
			graphinfo[name] = subs[2:]
f.close()

for pep in pepHash[0]:
	for info in pepHash[0][pep]:
		protein = info[pidx]
		if protein not in graphinfo:
			sys.exit("dbpep " + protein + " not found in graphinfo hash")
		ginfo = graphinfo[protein]
		edge = "\t".join(ginfo)
       	 	fw.write(pep+"\t"+protein+"\t"+str(info[eidx])+"\t"+des[0]+"\t"+str(info[cidx])+"\t"+edge+"\n")

#get protein sequence (for inference of peptide location in the protein from FGS)
protseq = {}
for pep in pepHash[1]:
	for info in pepHash[1][pep]:
		protein = info[pidx]
		protseq[protein] = ""
f = open(fgsfa, "r")
print "load", fgsfa
for aline in f:
	if aline[0] == '>':
		name = aline[1:-1]
	elif name in protseq:
		protseq[name] = aline[:-1]	
f.close()

#write fgs info to the output file
for pep in pepHash[1]:
	for info in pepHash[1][pep]:
		protein = info[pidx]
		seq = protseq[protein]
		a = len(pep)
		start = -1
		for p in range(len(seq)):
			if seq[p:p+a] == pep:
				start = p * 3 	
				break
		if start == -1:
			print "peptide", pep
			print "protein", protein
			print "seq",seq
			sys.exit("Warning: peptide not found in protein")
		start -= 1 #index from 0
		end = start + 3 * a
		subs = protein.split("_")
		edge = int(subs[7])
		if subs[-1] == '+':
			edge -= 1
       		fw.write(pep+"\t"+protein+"\t"+str(info[eidx])+"\t"+des[1]+"\t"+str(info[cidx])+"\t"+str(edge)+"\t"+str(start)+"\t"+str(edge)+"\t"+str(end)+"\t0"+"\n")

fw.close()
