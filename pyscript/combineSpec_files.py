import sys
from sys import argv
import re

if len(sys.argv) < 2:
	print sys.argv[0] + " at-least-one-spectral-file"
	sys.exit()

def spectraHash(inputs):
	titleHash = dict()
	specPep = dict()
	for input in inputs:
		print "process", input
		f = open(input, "r")
		for line in f:
			if line.startswith('#'):
                        	continue
	                lines = line.split('\t')
        	        title = lines[3]
			peptide = re.sub(r'[^a-zA-Z]', '', lines[9])
	                protein = lines[10]
	                if protein.startswith('XXX_') :
	                 	continue      
	                else:
				if (title in specPep) and (peptide not in specPep[title]):
					print "same spec", title, " different peptides", peptide, specPep[title]	
					titleHash[title] = titleHash[title] + 1
					tmp = specPep[title]
					tmp.append(peptide)
					specPep[title] = tmp
				else:
	                        	titleHash[title] = 1
					specPep[title] = [peptide,]
	f.close()
	return titleHash

mergedHash = spectraHash(argv[1:])
print "total identified spectra", len(mergedHash)
ambigious = 0
for key in mergedHash.keys():
	if mergedHash[key] > 1:
		ambigious += 1
print "ambigious", ambigious
