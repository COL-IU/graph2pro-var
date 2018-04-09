from sys import argv
filename = argv[1]
fgsfilename = argv[2]
outputfilename = argv[3]

queryHash = dict()
fgsHash = dict()
def parseFGS(fgsfilename):
	with open(fgsfilename) as f:
		for line in f:
			if line.startswith('>'):
				name = line.replace('>','')
				#seqID = name.rstrip()
				seqID = name.rstrip().split(' ')[0]
				#print seqID
				sequence = f.next().rstrip()
				fgsHash[seqID] = sequence
	f.close()
	print len(fgsHash)
	return(fgsHash)

fgsHash = parseFGS(fgsfilename)
resultHash = dict()
of = open(outputfilename,'w')

with open(filename) as f:
	for line in f:
		if line.startswith('#'):
			continue
		lines = line.split('\t')
		query = lines[0]
		mismatch = int(lines[4])
		identity = float(lines[2])
		if query in queryHash:
			continue
		else:
			queryHash[query] = 1
			#print query
			if query in fgsHash:
				#print query
				if (identity >= 70.0) and (identity <= 100.0):
					resultHash[query] = 1
					#print '>',lines[0],'\t',identity,'\t',mismatch,'\t',fgsHash[query]
					of.write('>'+query+'\n'+fgsHash[query]+'\n')
of.close()
#print len(resultHash)
