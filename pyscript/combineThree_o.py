#!/usr/bin/env python
from sys import argv

if len(argv) < 4:
        print('Usage: %s <mismatched> <unconcordant> <graph2pro>' % argv[0])
        exit(0)

mismatched = argv[1]
unconcondant = argv[2]
graph2pro = argv[3]

def spectraHash(input):
	titleHash = dict()
	with open(input) as f:
		for line in f:
			if line.startswith('#'):
                        	continue
	                lines = line.split('\t')
        	        title = lines[3]
	                protein = lines[10]
	                if protein.startswith('XXX_') :
	                 	continue      
	                else:
	                        titleHash[title] = 1
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
print len(mergedHash2)
#for eachID in mergedHash:
#	print eachID

