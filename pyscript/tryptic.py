#!/usr/bin/env python

import sys

infile, tryptic, nontryptic, cut = "", "", "", "no" 
minlen = 7

for idx in range(len(sys.argv)):
	if (sys.argv[idx] == "-i") and (len(sys.argv) > idx + 1):
		infile = sys.argv[idx + 1]
	elif (sys.argv[idx] == "-t") and (len(sys.argv) > idx + 1):
		tryptic = sys.argv[idx + 1]
	elif (sys.argv[idx] == "-n") and (len(sys.argv) > idx + 1):
		nontryptic = sys.argv[idx + 1]
	elif (sys.argv[idx] == "-c") and (len(sys.argv) > idx + 1):
		cut = sys.argv[idx + 1]
	elif (sys.argv[idx] == "-l") and (len(sys.argv) > idx + 1):
		minlen = sys.argv[idx + 1]
if not (infile and tryptic):
	sys.exit(sys.argv[0] + " -i input -t tryptic-pep-out <-n other-pep-out -c cut -l minlen>")

inf = open(infile, "r")
out1 = open(tryptic, "w")
if nontryptic:
	out2 = open(nontryptic, "w")

seqid = ""
for aline in inf:
	aline = aline.strip()
	if aline[0] == '>':
		seqid = aline 
	else:
		seqseq = aline
		slen = len(seqseq)
		leftmost = slen
		for p in range(slen):
			if (seqseq[p] == 'R') or (seqseq[p] == 'K'):
				leftmost = p
				break
		rightmost = 0
		for p in range(slen - 1, leftmost, -1):
			if (seqseq[p] == 'R') or (seqseq[p] == 'K'):
				rightmost = p
				break

		if rightmost - leftmost >= minlen:
			print >>out1, seqid
			if cut == "yes":
				print >>out1, seqseq[leftmost+1:rightmost+1]
			else:
				print >>out1, seqseq
		elif nontryptic:
			print >>out2, seqid 
			print >>out2, seqseq
inf.close()
out1.close()
if nontryptic:
	out2.close()
