#!/usr/bin/env python
#Nov 8 2018, Ye
import sys

def revcomplement(seq):
        code = {"A":"T", "T":"A", "C":"G", "G":"C", "a":"t", "t":"a", "c":"g", "g":"c", "X":"X", "N":"N"}
        comp = ""
        for base in seq[::-1]:
                comp += code[base]
        return comp

if len(sys.argv) < 4:
	sys.exit(sys.argv[0] + " contig-file old-output  new-output")
contigfile, oldinput, newoutput = sys.argv[1], sys.argv[2], sys.argv[3]

inf = open(oldinput, "r")
peplines = inf.readlines()
inf.close()

nucseq = {}
for aline in peplines:
	aline = aline.strip()
	subs = aline.split("\t")
	if subs[3] == 'FragGeneScan':
		pep, proid, pb, pe = subs[0], subs[1], int(subs[6]), int(subs[8])
		subs2 = proid.split("_")
		cid = "_".join(subs2[:-3])
		nucseq[cid] = ""

inf = open(contigfile, "r")
for aline in inf:
	aline = aline.strip()
	if aline[0] == '>':
		nucid = aline[1:]
	elif nucid in nucseq:
		nucseq[nucid] = aline
inf.close()

out = open(newoutput, "w")
codon={'TCA':'S','TCC':'S','TCG':'S','TCT':'S','TTC':'F','TTT':'F','TTA':'L','TTG':'L','TAC':'Y','TAT':'Y','TAA':'*','TAG':'*','TGC':'C','TGT':'C','TGA':'*','TGG':'W','CTA':'L','CTC':'L','CTG':'L','CTT':'L','CCA':'P','CCC':'P','CCG':'P','CCT':'P','CAC':'H','CAT':'H','CAA':'Q','CAG':'Q','CGA':'R','CGC':'R','CGG':'R','CGT':'R','ATA':'I','ATC':'I','ATT':'I','ATG':'M','ACA':'T','ACC':'T','ACG':'T','ACT':'T','AAC':'N','AAT':'N','AAA':'K','AAG':'K','AGC':'S','AGT':'S','AGA':'R','AGG':'R','GTA':'V','GTC':'V','GTG':'V','GTT':'V','GCA':'A','GCC':'A','GCG':'A','GCT':'A','GAC':'D','GAT':'D','GAA':'E','GAG':'E','GGA':'G','GGC':'G','GGG':'G','GGT':'G'};
for aline in peplines:
	aline = aline.strip()
	subs = aline.split("\t")
	if subs[3] != 'FragGeneScan':
		print >>out, aline
	else:
		pep, proid, pb, pe = subs[0], subs[1], int(subs[6]), int(subs[8])
		plen = len(pep)
		subs2 = proid.split("_")
		gb, ge, gdir = int(subs2[-3]), int(subs2[-2]), subs2[-1]
		l = int(subs2[3])
		cid = "_".join(subs2[:-3])
		nuc = nucseq[cid]
		nb = gb + pb + 1
		if gdir == '-':
			nb = l-ge+pb+1
			nuc = revcomplement(nuc)
		match = False
		nb_final = nb
		ptry = [0, 1, -1, 2, -2, 3, -3, 4, -4]
		#print "nb", nb
		for i in ptry: #try positions until find matched peptide 
			p = nb + i
			valid = False 
			mlen = 0
			for j in range(plen):
				k = p + j * 3;
				if len(nuc) < k+3: 
					break
				triplet = nuc[k:k+3]
				#print "triplet", triplet, "aa", codon[triplet]
				if triplet not in codon:
					break
				if codon[triplet] != pep[j]:
					break
				mlen += 1
			if mlen == plen:
				valid = True
				match = True
				nb_final = p
				#print "found match at try i", i
				break
		if not match:
			print "Warning: peptide not found in " + cid
		subs[6] = str(nb_final)
		subs[8] = str(nb_final + (pe - pb))
		print >> out, "\t".join(subs) 
out.close()
