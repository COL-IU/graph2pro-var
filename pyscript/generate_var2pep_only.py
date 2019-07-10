#!/usr/bin/env python
from sys import argv
import re

if len(argv) < 6:
        print('Usage: %s <mismatched> <unconcordant> <graph2pro> graph2pro_variant> <variant_only>' % argv[0])
        exit(0)

mismatched = argv[1]
unconcondant = argv[2]
graph2pro = argv[3]
graph2pro_variant = argv[4]
variant_only=argv[5]

def spectraHash(input):
	titleHash = dict()
	with open(input) as f:
		for line in f:
			if line.startswith('#'):
                        	continue
	                lines = line.split('\t')
        	        title = lines[3]
	                protein = lines[10]
			peptide =re.sub(r'[^a-zA-Z]', '', lines[9])
	                if protein.startswith('XXX_') :
	                 	continue      
	                else:
	                        titleHash[peptide] = 1
	f.close()
	return titleHash

def merge_two_dicts(x, y):
    '''Given two dicts, merge them into a new dict as a shallow copy.'''
    z = x.copy()
    z.update(y)
    return z

titleHash1 = spectraHash(mismatched)
titleHash2 = spectraHash(unconcondant)
mergedHash = merge_two_dicts(titleHash1,titleHash2)
titleHash3 = spectraHash(graph2pro)
mergedHash2 = merge_two_dicts(mergedHash,titleHash3)

gv = open(graph2pro_variant,'w')
#print len(mergedHash2)
for eachID in mergedHash2:
	gv.write(eachID+'\n')
gv.close()

vo = open(variant_only,'w')
for eachID in mergedHash :
	if eachID not in titleHash3 :
		vo.write(eachID+'\n')
vo.close()
