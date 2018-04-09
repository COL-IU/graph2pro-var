import sys

if len(sys.argv) < 3:
	sys.exit("need inputfile outfile")

input_file = sys.argv[1]
output_file = sys.argv[2]

ow = open(output_file,'w')
inf = open(input_file, 'r')
valid = 1
for aline in inf:
	if aline[0] == '>':
		subs = aline[1:-1].split()
		seqid = subs[0]
		if "'" not in seqid :
			ID = seqid.split(':')[0].replace(';','')
			ow.write('>'+ID+'\n')
			valid = 1
		else:
			valid = 0
	elif valid == 1:
		ow.write(aline)
ow.close()
inf.close()
