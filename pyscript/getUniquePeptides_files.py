#!/usr/bin/env python
from sys import argv
import sys
import re

def getUniqueHash(filenames, pephash0={}, spechash0={}):
	pepHash = dict()
	specHash = dict()
	scoreHash = dict()
	for afile in filenames:
		f = open(afile, "r")
		f.next()
		for line in f:
			lines = line.split('\t')
			peptide = re.sub(r'[^a-zA-Z]', '', lines[9])
			proteins = lines[10].split('(')[0]
			title = lines[3]
			score = int(lines[12])
			if proteins.startswith('XXX'):
				continue;
			#if multiple peptides were identified from a single spectral, only one is kept
			if (title in spechash0) or (title in specHash): #already identified spectrum
				if (title in specHash) and (score > scoreHash[title]): #keep the one with higher score
					specHash[title] = peptide
					scoreHash[title] = score
			else:
				if peptide not in pephash0:
					if peptide in pepHash:
						pepHash[peptide]=pepHash[peptide]+1
					else:
						pepHash[peptide] = 1
				specHash[title] = peptide	
				scoreHash[title] = score
		f.close()
	return pepHash, specHash
def writeTofile(hashinput,wfile,tag):
	for i in hashinput:
		key = i
		value = hashinput[key]
		wfile.write(key+'\t'+str(value)+'\t'+tag+'\n')

infiles, varfiles = [], []
outputfile = ""
brief=False
for idx in range(len(argv)):
	if (argv[idx] == '-i') and len(argv) > idx + 1:
		infiles.append(argv[idx+1])	
	elif (argv[idx] == '-v') and len(argv) > idx + 1:
		varfiles.append(argv[idx+1])	
	elif (argv[idx] == '-o') and len(argv) > idx + 1:
		outputfile = argv[idx + 1]	
	elif (argv[idx] == '-b'):
		brief = True
if len(infiles) < 1:
	print sys.argv[0] + " -i in-file (at least one) [-v spec-file] [-o ouptut]"
	print "an output file name, and at least one tsv file need to be provided"
	sys.exit()

pephash, spechash = getUniqueHash(infiles)
if outputfile:
	out = open(outputfile, "w")
	writeTofile(pephash,out,"graph2pro")
numbers = []
if brief:
	numbers.extend([str(len(spechash)), str(len(pephash))])
else:
	print "+".join(infiles), "spec", len(spechash), "unique-peptide", len(pephash)
if len(varfiles) > 0:
	pephash_only, spechash_only = getUniqueHash(varfiles, pephash, spechash)
	if brief:
		numbers.extend([str(len(spechash_only)), str(len(pephash_only))])
	else:
		print "+".join(varfiles) + "-only", "spec", len(spechash_only), "unique-peptide", len(pephash_only)
	allfiles = infiles[:]
	allfiles.extend(varfiles)
	if brief:
		numbers.extend([str(len(spechash) + len(spechash_only)), str(len(pephash) + len(pephash_only))])
	else:
		print "+".join(allfiles) + "-all", "spec", len(spechash) + len(spechash_only), "unique-peptide", len(pephash) + len(pephash_only)
	if outputfile:
		writeTofile(pephash_only, out,"var2pep")
if outputfile:
	out.close()
if brief:
	print " ".join(numbers)
