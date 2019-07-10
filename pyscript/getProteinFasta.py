#!/usr/bin/env python
from sys import argv
from Bio import SeqIO

if len(argv) < 2:
        print('Usage: %s <protein> <input_fastafile>' % argv[0])
        exit(0)

def readFasta(fasta):
        handle = open(fasta, "rU")
        record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        handle.close()
        return record_dict

protein = argv[1]
fasta = argv[2]
fastaHash = readFasta(fasta)

with open(protein) as pf:
        for line in pf:
                lines = line.split('\t')
                protein = lines[0]
                if protein in fastaHash:
                        print '>'+protein+'\n'+fastaHash[protein].seq

