from sys import argv
from Bio import SeqIO
fastafile = argv[1]


for seq_record in SeqIO.parse(fastafile,"fasta"):
	#seqID = seq_record.id
	seqDes = seq_record.description
	seqDes_2 = '_'.join(seqDes.split(' '))
	seqString = str(seq_record.seq)
	print '>'+seqDes_2+'\n'+seqString
