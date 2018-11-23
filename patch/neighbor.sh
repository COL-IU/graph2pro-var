awk 'FNR==NR{ if($4=="FragGeneScan") { n=split($2, t, "_"); cid=t[8]; cidlist[cid] = $1; pidlist[$2] = 1; } next; } {s=substr($1,1,1); if (s==">") { m=split($1,p,"_"); pid=substr($1,1); cid=p[8]; if((cid in cidlist) && (pid in pidlist==0)) { valid=1; print $0; } else { valid=0; }} else if(valid) print $0; }' FS="\t" Bst45.step1.output Bst45.contig.fgs.fasta > Bst45.neighbor.fasta

java -Xmx16g -jar /u/yye/Metaproteomics/MSGF+/MSGFPlus.jar -s ../MS/2016_Jan_12_QE2_45.mgf -o Bst45.neighbor.mzid -d Bst45.neighbor.fasta -inst 1 -t 15ppm -ti -1,2 -mod /u/yye/Metaproteomics/MSGF+/Mods_normal.txt -ntt 1 -tda 1 -maxCharge 7 -minCharge 1 -addFeatures 1 -n 1

java -Xmx16g -cp /u/yye/Metaproteomics/MSGF+/MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i Bst45.neighbor.mzid  -showDecoy 1
