from sys import argv
import re

def getUniqueHash(filename):
	pepHash = dict()
	proHash = dict()
	with open(filename) as f:
		f.next()
		for line in f:
			lines = line.split('\t')
			peptide = re.sub(r'[^a-zA-Z]', '', lines[9])
			proteins = lines[10].split(';')
			for eachPro in proteins:
				eachPurePro = eachPro.split('(')[0]
				if eachPurePro.startswith('XXX'):
					continue;
				if eachPurePro in pepHash:
					pepHash[eachPurePro] = pepHash[eachPurePro]+1
				else:
					pepHash[eachPurePro] = 1

	return(pepHash)
def writeTofile(hashinput,wf):
	wfile = open(wf,'w')
	for i in hashinput:
		key = i
		value = hashinput[key]
		wfile.write(key+'\t'+str(value)+'\n')
	wfile.close()

hash1 = getUniqueHash(argv[1])
print len(hash1)
writeTofile(hash1,argv[1]+'.uniquePro')
#hash2 = getUniqueHash(argv[2])
#hash3 = getUniqueHash(argv[3])
#hash4 = getUniqueHash(argv[4])

#intersection = set(hash1.keys()) & set(hash2.keys())

#print len(hash1),len(hash2),len(intersection)
