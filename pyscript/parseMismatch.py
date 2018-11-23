#!/usr/bin/env python
#Ye Nov 7 2018

import sys

if(len(sys.argv) < 4):
	sys.exit(sys.argv[0] + " similarty-research-file fgs-file outputfile <identity cutoff>")

filename = sys.argv[1]
fgsfilename = sys.argv[2]
outputfilename = sys.argv[3]
if(len(sys.argv) > 4):
	identcut = float(sys.argv[4])
	if identcut < 1.0:
		identcut *= 100
else:
	identcut = 70

def processOne(lines):
	maxmatchlen = 0
	for aline in lines:
		#print aline
		subs = aline.split('\t')
		identity = float(subs[2])
		matchlen = int(subs[3])
		#if (identity == 100.0) and (matchlen >= 30):
		#if (identity >= 95.0) and (matchlen >= 30):
		#	maxmatchlen = 0
		#	break #already presented in the identified proteins
		if (identity >= identcut) and (matchlen >= 20): #add matchlen cutoff
			if matchlen > maxmatchlen:
				maxmatchlen = matchlen
	#print "maxmatchlen: ", maxmatchlen 
	#raw_input("type enter to continue...")
	return maxmatchlen 

queryHash = dict()
lastlines = []
lastquery = ""
with open(filename) as f:
	for line in f:
		if line.startswith('#'):
			continue
		subs = line.split('\t')
		if subs[0] != lastquery:
			if lastlines:
				mlen = processOne(lastlines)
				if mlen > 0:
					queryHash[lastquery] = mlen
			lastquery = subs[0]
			lastlines = []
		lastlines.append(line)
	if lastlines:
		mlen = processOne(lastlines)
		if mlen > 0:
			queryHash[lastquery] = mlen
	f.close()
print "queryHash done"

of = open(outputfilename,'w')
matchlen = 0
tryptic = {}
minlen = 7
with open(fgsfilename) as f:
	for line in f:
		line = line.strip()
		if line[0] == '>':
			name = line[1:]
			matchlen = 0
			if name in queryHash: 
				matchlen = queryHash[name]
		elif matchlen > 0:
			seqseq = line
			slen = len(seqseq)
                	leftmost = slen
			#only outputs those that potentially contain a tryptic peptide
	                for p in range(slen):
				if (seqseq[p] == 'R') or (seqseq[p] == 'K'):
       	                        	leftmost = p
                                	break
                	rightmost = 0
                	for p in range(slen - 1, leftmost, -1):
                        	if (seqseq[p] == 'R') or (seqseq[p] == 'K'):
                                	rightmost = p
                                	break
                	if rightmost - leftmost < minlen:
				continue
			if seqseq[leftmost:rightmost] not in tryptic:
				of.write('>'+name+'\n'+line+'\n')
				tryptic[seqseq[leftmost:rightmost]] = name
	f.close()
of.close()
