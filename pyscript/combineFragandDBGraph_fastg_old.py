import sys
import re
from sys import argv
#from utils import *
#from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna

def revcomplement(seq):
        code = {"A":"T", "T":"A", "C":"G", "G":"C", "a":"t", "t":"a", "c":"g", "g":"c", "X":"X", "N":"N"}
        comp = ""
        for base in seq[::-1]:
                comp += code[base]
        return comp

totalLen = len(argv)
if totalLen != 8:
	print ("The argv should be <FGS> <FGSNA> <FGSPro> <DBGraph> <DBGraphDatabase> <contig> <OutputFile>")
	sys.exit()

FragGeneScan = argv[1]
FragGeneScanNA = argv[2]
FragGeneScanPro = argv[3]
DBGraph = argv[4]
DBGraphDatabase = argv[5]
contig = argv[6]
CombinedFile = argv[7]

def readDBGraphHead(DBGraphDatabase):
	head = dict()
	with open(DBGraphDatabase) as f:
		for line in f:
			if line.startswith('>'):
				line = line.rstrip()
				line = line.replace('>','')
				lines = line.split(' ')
				seq = f.next().rstrip()
				head[lines[0]] = line+" "+seq
	f.close()
	return head
def readFragGeneScan(FragGeneScanDatabase) :
	fastaHash = dict()
	with open(FragGeneScanDatabase) as f:
		for line in f:
			if line.startswith('>'):
				head = line.rstrip().split(' ')[0].replace('>','')
				#print head
				seq = f.next()
				seq = seq.rstrip()
				fastaHash[head] = seq
	f.close()
	return fastaHash

def readFastgContig(contig_fastg) :
        fastaHash = dict()
        with open(contig_fastg) as f:
                for line in f:
                        if line.startswith('>'):
                                head = '_'.join(line.replace('>','').split('_')[0:1])
                                #print head
                                seq = f.next()
                                seq = seq.rstrip()
                                fastaHash[head] = seq
        f.close()
        return fastaHash

class pepStruc(object):
    peptide = None
    bestEvalue = None
    peptideCounts = None
    protein = None
    edge1 = None
    start = None
    edge2 = None
    end = None
    stopCodon = None
    source = None
    strain = '+'

DBGraphHead = readDBGraphHead(DBGraphDatabase)
FGSNAFasta = readFragGeneScan(FragGeneScanNA)
FGSProFasta = readFragGeneScan(FragGeneScanPro)
contigFasta = readFastgContig(contig)

fw = open(CombinedFile,'w')
fw.write('#peptide\tprotein\tbestEvalue\tsource\tpeptideCounts\tedge1\tstart\tedge2\tend\tstopCodon\n')
peptideHash = dict()
with open(FragGeneScan) as f:
	for line in f:
		line = line.rstrip()
		if line.startswith("#"):
			continue
		else:
			lines = line.split('\t')
			peptide = re.sub('[^a-zA-Z]', '', lines[9])
			
			if peptide in peptideHash :
                                pepS = peptideHash[peptide]
                                pepS.peptideCounts += 1
                                if evalue < pepS.bestEvalue :
                                        pepS.bestEvalue = evalue
                                peptideHash[peptide] = pepS
				continue
			protein = lines[10].split('(')[0]
			#print protein
			if protein.startswith('XXX_'):
                                continue
			#print peptide,'\t',protein
			geneSequence = FGSNAFasta[protein][0:100]

			proteins = protein.split('_')
			leng = len(proteins)
			geneid = "_".join(proteins[0:1])			

			#for fastg input
			edge1 = int(proteins[leng-4])
			print edge1
                        start = int(proteins[leng-3])
                        edge2 = edge1
                        end = int(proteins[leng-2])
			strand = proteins[leng-1]

			#strand = proteins[3]
		        #edge1 = int(proteins[0])*2 - 1
                        #start = int(proteins[1])*2 - 1
                        #edge2 = int(proteins[0])*2 - 1
                        #end = int(proteins[2])*2 -1 
			
			#print geneSequence
			contigSequence = contigFasta[geneid]
			geneLocation = contigSequence.find(geneSequence)
			if geneLocation < 0 :
				#my_dna = Seq(contigSequence,generic_dna)
				#my_dna_rc = my_dna.reverse_complement()
				my_dna_rc = revcomplement(contigSequence)
				geneLocation = my_dna_rc.find(geneSequence)
			#geneLocation = findGeneLocation(contig,proteins[0],geneSequence)
			#print geneLocation
			if '+' in strand :
				edge1 = edge1 - 1
				edge2 = edge2 - 1
                        stopCodon = 0
                        evalue = lines[14]
			
			#if protein in FGSFasta :
			FGSseq = FGSProFasta[protein]
			location1 = FGSseq.find(peptide)*3+geneLocation
			location2 = location1+len(peptide)*3
			#print location1,location2
			start = location1
			end = location2
			if 'post=-' in protein :
				stopCodon = 1
			tmp = pepStruc()
			tmp.peptide = peptide
			tmp.protein = protein
			tmp.bestEvalue = evalue
			tmp.peptideCounts = 1
			tmp.edge1 = edge1
			tmp.start = start
			tmp.edge2 = edge2
			tmp.end = end
			tmp.stopCodon = stopCodon
			tmp.source = 'FragGeneScan'
			tmp.strand = strand
			print peptide,protein,edge1,start,edge2,end,stopCodon
			peptideHash[peptide] = tmp
f.close()
#for eachPep in peptideHash :
#	pepS = peptideHash[eachPep]
#	fw.write(pepS.peptide+"\t"+pepS.protein+"\t"+pepS.bestEvalue+"\t"+pepS.source+"\t"+str(pepS.peptideCounts)+"\t"+str(pepS.edge1)+"\t"+str(pepS.start)+"\t"+str(pepS.edge2)+"\t"+str(pepS.end)+"\t"+str(pepS.stopCodon)+"\n")


peptideHash2 = dict()
with open(DBGraph) as f:
	for line in f:
		line = line.rstrip()
		if line.startswith('#'):
			continue
		else:
			lines = line.split('\t')
			peptide = re.sub(r'[^a-zA-Z]', '', lines[9])
			print peptide
			tmpL = lines[10]
			evalue = lines[13]
			proteins = tmpL.split('(')
			protein = proteins[0]
			#print protein
			if protein.startswith('XXX_'):
				continue
			proteinDes = DBGraphHead[protein]
			#print peptide,'\t',protein,'\t',proteinDes
			#if protein in DBGraphHead :
			#	proteinDes = DBGraphHead[protein]
			#print proteinDes
			proteinDes_s = proteinDes.split(' ');
			edge1 = int(proteinDes_s[2])
			start = int(proteinDes_s[3])
			edge2 = int(proteinDes_s[4])
			end = int(proteinDes_s[5])
			stopCodon = int(proteinDes_s[6])
			orig_pep = peptide
			if orig_pep in peptideHash2 :
				print orig_pep
                                pepS = peptideHash2[orig_pep]
				if edge1 != pepS.edge1 :
					print peptide,protein,edge1,start,edge2,end
					print pepS.peptide,pepS.protein,pepS.edge1,pepS.start,pepS.edge2,pepS.end
                                pepS.peptideCounts += 1
                                if evalue < pepS.bestEvalue :
                                        pepS.bestEvalue = evalue
                                peptideHash2[peptide] = pepS
                        else:
				print orig_pep
                                tmp = pepStruc()
                                tmp.peptide = orig_pep
                                tmp.protein = protein
                                tmp.bestEvalue = evalue
                                tmp.peptideCounts = 1
                                tmp.edge1 = edge1
                                tmp.start = start
                                tmp.edge2 = edge2
                                tmp.end = end
                                tmp.stopCodon = stopCodon
                                tmp.source = 'DBGraph2Pro'
				tmp.strand='+'
                                peptideHash2[peptide] = tmp

print len(peptideHash2)
for eachPep in peptideHash2 :
        pepS = peptideHash2[eachPep]
	fw.write(pepS.peptide+"\t"+pepS.protein+"\t"+pepS.bestEvalue+"\t"+pepS.source+"\t"+str(pepS.peptideCounts)+"\t"+str(pepS.edge1)+"\t"+str(pepS.start)+"\t"+str(pepS.edge2)+"\t"+str(pepS.end)+"\t"+str(pepS.stopCodon)+"\n")

for eachPep in peptideHash :
        pepS = peptideHash[eachPep]
	if eachPep in peptideHash2 :
		pepS2 = peptideHash2[eachPep]
		edgeFGS = pepS.edge1
		FGSstart = pepS.start
		FGSend = pepS.end
		edgeArray = [pepS2.edge1,pepS2.edge2]
		coords = [pepS2.start,pepS2.end]
		if pepS2.edge1 != pepS2.edge2 :
			if edgeFGS in edgeArray:
				if (FGSstart in coords) or (FGSend in coords) :
					continue 
        	fw.write(pepS.peptide+"\t"+pepS.protein+"\t"+pepS.bestEvalue+"\t"+pepS.source+"\t"+str(pepS.peptideCounts)+"\t"+str(pepS.edge1)+"\t"+str(pepS.start)+"\t"+str(pepS.edge2)+"\t"+str(pepS.end)+"\t"+str(pepS.stopCodon)+"\n")
	else:
		fw.write(pepS.peptide+"\t"+pepS.protein+"\t"+pepS.bestEvalue+"\t"+pepS.source+"\t"+str(pepS.peptideCounts)+"\t"+str(pepS.edge1)+"\t"+str(pepS.start)+"\t"+str(pepS.edge2)+"\t"+str(pepS.end)+"\t"+str(pepS.stopCodon)+"\n")
fw.close()

#for eachPep in peptideHash:
#	if eachPep in peptideHash2 :
#		dbgraph = peptideHash2[eachPep]
#                pep = peptideHash[eachPep]
		#print 'DBGraph',dbgraph.peptide,dbgraph.edge1,dbgraph.start,dbgraph.edge2,dbgraph.end,dbgraph.stopCodon,dbgraph.strand
                #print 'FragGeneScan:',pep.peptide,pep.edge1,pep.start,pep.edge2,pep.end,pep.stopCodon,pep.strand
