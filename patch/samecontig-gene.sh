#awk 'FNR==NR{ n=split($1, t, "_"); cid=t[8]; cidlist[cid] = $1; next; } {s=substr($1,1,1); if (s==">") { m=split($1,p,"_"); cid=p[8]; if(cid in cidlist) { valid=1; print $0; } else { valid=0; }} else if(valid) print $0; }' FS="\t" fgs-id.list Bst45.contig.fgs.faa > samecontig-gene.faa

awk 'FNR==NR{ if($4=='FragGeneScan') { n=split($2, t, "_"); cid=t[8]; cidlist[cid] = $1; pidlist[$2] = 1; } next; } {s=substr($1,1,1); if (s==">") { m=split($1,p,"_"); pid=substr($1,1); cid=p[8]; if(cid in cidlist) { if(pid in pidlist) { valid = 0; } else {valid=1; print $0; }} else { valid=0; }} else if(valid) print $0; }' FS="\t" Bst45.step1.output Bst45.contig.fgs.faa > Bst45.oneneib.faa

