awk 'FNR==NR{if($4 == "FragGeneScan"){  n=split($2, t, "_"); id=t[8]; idlist[id] = $2; } next; } {s=substr($1,1,1); if (s==">") { m=split($1,p,"_"); id=p[8]; if(id in idlist) { valid=1; print $0; } else { valid=0; }} else if(valid) print $0; }' FS="\t" Bst45.step1.output Bst45.contig.fasta > Bst45.contig.sel

python patchpos.py Bst45.contig.sel Bst45.step1.output Bst45.step1.output.fixed

mv Bst45.step1.output.fixed Bst45.step1.output

awk 'FNR==NR{idlist[$2]; next} {s=substr($1,1,1); if (s==">") { id=substr($1, 2); if (id in idlist) { valid = 1; print ">" id; } else { valid = 0; }} else if (valid) print $0; }' FS="\t" Bst45.step1.output Bst45.contig.fgs.fasta > Bst45.step1.fasta.1

/u/yye/Metaproteomics/Graph2Pro/DBGraphPep2Pro -s /home/team/metaproteomics/ocean/Assembly-megahit/Bst/k99.contigs.fastg -f -p /home/team/metaproteomics/ocean/Bst45/Bst45.step1.output -o Bst45.step1.fasta.2 -k 99 -d 5

cat Bst45.step1.fasta.1  Bst45.step1.fasta.2 > Bst45.step1.fasta.12

cd-hit -i Bst45.step1.fasta.12 -c 1.0 -o Bst45.step1.fasta -d 0
