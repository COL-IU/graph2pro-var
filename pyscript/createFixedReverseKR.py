import sys
from sys import argv

def revfix(seq):
	seq_rev = ""
	if seq[len(seq)-1] == 'K' or seq[len(seq) -1] == 'R' :
		seq_2 = seq[0:len(seq)-1]
		seq_rev = seq_2[::-1]+seq[len(seq)-1]
	else:
		seq_rev = seq[::-1]
	return seq_rev

def processone(f, des, seq):
	f.write(">" + des + "\n" + seq + "\n")
	seq_rev = revfix(seq)
	f.write(">XXX_" + des + "\n" + seq_rev + "\n")

if len(sys.argv) < 2:
	sys.exit("need inputfile")

filename = argv[1]
write_filename = argv[1].replace('.fasta','.fixedKR.fasta')
print write_filename
out = open(write_filename,'w')
inf = open(argv[1], "r")
des, seq = "", ""
for aline in inf:
	if aline[0] == '>':
		if des:
			processone(out, des, seq)
		des = aline[1:-1]
		seq = ""
	else:	
		seq += aline.strip()
if seq:
	processone(out, des, seq)	
out.close()
inf.close()

