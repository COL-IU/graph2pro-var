#!/bin/bash

#make sure perl, python, java are available when this pipeline runs

#don't change the following
#$0 only works for direct submission
script=$(basename $(readlink -nf $0))
pgm_d=$(dirname $(readlink -nf $0))
#if use qsub, the folder info needs to be passed along using the parameter 
if [ ! -z $pgmdir ]; then
   pgm_d=$pgmdir 
fi

parfile=$1
if [ ! -z $par ]; then
  parfile=$par
fi
if [ -z $parfile ]; then
  echo "Error: parameter file not given";
  exit;
fi
if [ -f $parfile ]; then
   echo "$parfile found -- nice"
else
   echo "$parfile NOT found -- pipeline terminated"
   exit
fi
#get parameters from input parfile
IFS=$'\n'; set -f; par=($(<$parfile))
n=${#par[@]}
if [ $n -lt 5 ]; then
   echo "Error: wrong parameter file"
   exit
fi
thread=8
ram="32g"
fdr=0.01
ide=70
cascaded="no"
for apar in "${par[@]}"; do
   IFS='=' read -r -a item <<< "$apar" 
   case "${item[0]}" in
     id) exp_d=${item[1]};;
     kmer) kmer=${item[1]};;
     fastg) gnm=${item[1]};;
     ms) mgf=${item[1]};;
     reads) read_f=${item[1]};;
     thread) thread=${item[1]};;
     memory) ram=${item[1]};;
     fdr) fdr=${item[1]};;
     ide) ide=${item[1]};;
     cascaded) cascaded=${item[1]};;
   esac
done

if [ -f $gnm ]; then
   echo "$gnm found"
else
   echo "$gnm not found"
   exit
fi
if [ -f $mgf ]; then
   echo "$mgf found"
else
   echo "$mgf not found"
   exit
fi

if [ -z $read_f ]; then
  echo "Since NO reads files are given, only graph2pro will be run (var2pep step will be skipped)"
fi

#programs
dbgraph2pro="$pgm_d/Graph2Pro/DBGraph2Pro"
dbpep2pro="$pgm_d/Graph2Pro/DBGraphPep2Pro"
fgs="$pgm_d/FragGeneScan1.30/run_FragGeneScan.pl"
rap_d="$pgm_d/RAPSearch/bin"
rap="$rap_d/rapsearch"
msgf="$pgm_d/MSGF+/MSGFPlus.jar"
bowtie_d="$pgm_d/bowtie2-2.3.3.1"
programs=($dbgraph2pro $fgs $rap $msgf "$bowtie_d/bowtie2")
echo "check programs..."
for prog in ${programs[@]}; do
   if [ -f ${prog} ]; then
      echo "  $prog found -- nice"
   else
      echo "  $prog not found -- the pipeline is terminated"
      exit
   fi
done


#prepare Graph2Pep
if [ -f $exp_d.graph2pep.fasta ]; then
   echo "Graph2pep already done -- skip this step"
else
   echo "Now run Graph2Pep..."
   $dbgraph2pro -s $gnm -S -k $kmer -o $exp_d.graph2pep.fasta
   python $pgm_d/pyscript/createFixedReverseKR.py $exp_d.graph2pep.fasta
fi

#prepare contig file and fgs
if [ -f $exp_d.contig.fasta ]; then
	echo "$exp_d.contig.fasta found"
else
	echo "Now prepare $exp_d.contig.fasta file..."
	python $pgm_d/pyscript/fastg2fasta.py $gnm $exp_d.contig.fasta 
fi

if [ -f $exp_d.contig.fgs.fasta ]; then
   echo "FGS prediction for contigs already done -- skip this step"
else
   echo "Now run FGS for contigs..."
   $fgs -genome=$exp_d.contig.fasta -out=$exp_d.contig.fgs -complete=0 -train=complete -thread=$thread
   mv $exp_d.contig.fgs.faa $exp_d.contig.fgs.fasta
fi

#MSGF+ search against contigs
#first search
if [ -s $exp_d.fgs.mzid ]; then
    echo "MSGF+ against contigs already done -- skip this step"
else
    rm -f *.contig.fgs.revCat.*
    echo "Now run MSGF+ against contigs..."
    echo "java -Xmx$ram -jar $msgf -s $mgf -o $exp_d.fgs.mzid -d $exp_d.contig.fgs.fasta -inst 1 -t 15ppm -ti -1,2 -mod $pgm_d/MSGF+/Mods_normal.txt -ntt 1 -tda 1 -maxCharge 7 -minCharge 1 -addFeatures 1 -n 1 -thread $thread"
    java -Xmx$ram -jar $msgf -s $mgf -o $exp_d.fgs.mzid -d $exp_d.contig.fgs.fasta -inst 1 -t 15ppm -ti -1,2 -mod $pgm_d/MSGF+/Mods_normal.txt -ntt 1 -tda 1 -maxCharge 7 -minCharge 1 -addFeatures 1 -n 1 -thread $thread
    java -Xmx$ram -cp $msgf edu.ucsd.msjava.ui.MzIDToTsv -i $exp_d.fgs.mzid  -showDecoy 1
fi

python $pgm_d/pyscript/parseFDR_o.py $exp_d.fgs.tsv $fdr


#MSGF+ search against graph2pep
if [ -s $exp_d.graph2pep.mzid ]; then
    echo "MSGF+ against graph2pep already done -- skip this step"
else
    rm -r -f *.graph2pep.mzid
    echo "Now run MSGF+ against graph2pep..."
    java -Xmx$ram -jar $msgf -s $mgf -o $exp_d.graph2pep.mzid -d $exp_d.graph2pep.fixedKR.fasta -inst 1 -t 15ppm -ti -1,2 -mod $pgm_d/MSGF+/Mods_normal.txt -ntt 1 -tda 0 -maxCharge 7 -minCharge 1 -addFeatures 1 -n 1 -thread $thread
    java -Xmx$ram -cp $msgf edu.ucsd.msjava.ui.MzIDToTsv -i $exp_d.graph2pep.mzid -showDecoy 1
fi

#run DBGraphPep2Pro and prepare *.step1.fasta
if [ -s $exp_d.step1.output ]; then
    echo "combineFragDBGraph already done -- skip this step"
else
    rm -f *tsv.0.10.tsv
    echo "Now run combineFragDBGraph..."
    python $pgm_d/pyscript/parseFDR_o.py $exp_d.fgs.tsv 0.10
    python $pgm_d/pyscript/parseFDR_o.py $exp_d.graph2pep.tsv 0.10
    #prepare step1.output 
    python $pgm_d/pyscript/step1output.py $exp_d.fgs.tsv.0.10.tsv $exp_d.contig.fgs.fasta $exp_d.graph2pep.tsv.0.10.tsv $exp_d.graph2pep.fasta $exp_d.step1.output.old
    #patch positions
    python $pgm_d/pyscript/patchpos.py $exp_d.contig.fasta $exp_d.step1.output.old $exp_d.step1.output
    rm -f $exp_d.step1.output.old
    #end patch
fi

#prepare Graph2Pro database
if [ -s $exp_d.graph2pro.fasta ]; then
   echo "GraphPep2Pro already done -- skip this step"
else
   rm -f $exp_d.step1.fasta $exp_d.graph2pro.mzid $exp_d.uncordant.rapsearch.fasta $exp_d.mismatched.mzid $exp_d.uncordant.mzid 
   echo "Now run GraphPep2Pro..."
   echo "$dbpep2pro -s $gnm -f -p $exp_d.step1.output -o $exp_d.step1.fasta -k $kmer -d 5"
   $dbpep2pro -s $gnm -f -p $exp_d.step1.output -o $exp_d.step1.fasta.2 -k $kmer -d 5
   #merge proteins from contigs & proteins from graph and remove redundancy (using cd-hit)
   awk 'FNR==NR{idlist[$2]; next} {s=substr($1,1,1); if (s==">") { id=substr($1, 2); if (id in idlist) { valid = 1; print ">" id; } else { valid = 0; }} else if (valid) print $0; }' FS="\t" $exp_d.step1.output $exp_d.contig.fgs.fasta > $exp_d.step1.fasta.1
   cat $exp_d.step1.fasta.1  $exp_d.step1.fasta.2 > $exp_d.step1.fasta.12
   $pgm_d/cd-hit/cd-hit -i $exp_d.step1.fasta.12 -c 1.0 -o $exp_d.step1.fasta -d 0
   mv $exp_d.step1.fasta $exp_d.graph2pro.fasta
fi

#search against graph2pro
if [ -s $exp_d.graph2pro.mzid ]; then
   echo "Search against graph2pro db done -- skip this step"
else
   rm -f $exp_d.graph2pro.mzid
   rm -f $exp_d.graph2pro.revCat*
   echo "Now search against graph2pro db..."
   java -Xmx$ram -jar $msgf -s $mgf -o $exp_d.graph2pro.mzid -d $exp_d.graph2pro.fasta -inst 1 -t 15ppm -ti -1,2 -mod $pgm_d/MSGF+/Mods_normal.txt -ntt 1 -tda 1 -maxCharge 7 -minCharge 1 -addFeatures 1 -n 1 -thread $thread
   java -Xmx$ram -cp $msgf edu.ucsd.msjava.ui.MzIDToTsv -i $exp_d.graph2pro.mzid -showDecoy 1
   python $pgm_d/pyscript/parseFDR_o.py $exp_d.graph2pro.tsv $fdr
   python $pgm_d/pyscript/parseFDR_o_peptide.py $exp_d.graph2pro.tsv $fdr
fi

#get rejected spectra after graph2pro search, if "cascaded" option is turned on
if [ "${cascaded}" = "yes" ]; then
	#only exclude spectra matching proteins in the Graph2Pro database, but not the reverse (!=XXX)
	#use ~ /END IONS/ (== "END IONS" don't work with mgf files ending with ^M
	awk 'FNR==NR{protidstart=substr($11, 1, 3); if(protidstart!="XXX") idlist[$4]; next} {if($0 ~ /END IONS/) { if(id in idlist == 0) { for(x=1;x<=add;x++) print lines[x]; print "END IONS";} } if($0 ~ /BEGIN IONS/) { lines[1] = $0; add=1; } else { n=split($0, str, "="); if(str[1]=="TITLE") id=str[2]; add = add + 1; lines[add]=$0; }}' $exp_d.graph2pro.tsv.${fdr}.tsv ${mgf} > rejected_g2pX.mgf
	mgf="rejected_g2pX.mgf"
fi

if [ -z $read_f ]; then
   exit
fi

#var2pep
echo "Now run steps to collect read_f containing variants..."

if [ -f $exp_d.uncordant.fasta ]; then
   echo "Unmapped reads already extracted -- skip this step"
else
   #mapping
   if [ -f $exp_d.contig.bw2.1.bt2 ]; then
      echo "Bowtie2 database found -- skip this step"
   else
      echo "Now prepare bowtie2 database..."
      $bowtie_d/bowtie2-build $exp_d.contig.fasta $exp_d.contig.bw2
   fi

   if [ -s $exp_d.sam ]; then
      echo "Reads mapping to contigs already done -- skip this step"
   else
      echo "Now run bowtie2..."
      cmd="$bowtie_d/bowtie2 -x $exp_d.contig.bw2 ${read_f} --no-head --no-unal --un-conc $exp_d.un-conc.fq --un $exp_d.un.fq -S $exp_d.sam -p $thread"
      echo $cmd
      #note: use eval $cmd (direct call doesn't work!!)
      eval $cmd
      if [ -s $exp_d.sam ]; then
         echo "$exp_d.sam looks good"
      else
         echo "empty $exp_d.sam -- the pipeline is terminated"
         exit
      fi
   fi

   #extract reads
   awk -F'\t' '{if(($17 != "NM:i:0") && ($17 != "NM:i:1")) print ">" $1 " " $2 " " $6 " " $17 " " $18 "\n" $10}' $exp_d.sam > $exp_d.mismatched.higherthan1.fasta
   cat ${exp_d}.un-conc.1.fq ${exp_d}.un-conc.2.fq ${exp_d}.un.fq > ${exp_d}.uncordant.fq
   awk 'NR % 4 == 1 {print ">" substr($0,2) } NR % 4 == 2 {print $0}' ${exp_d}.uncordant.fq > ${exp_d}.uncordant.fasta
fi

if [ -f $exp_d.uncordant.fgs.ffn ]; then 
   echo "FGS of variant reads already done -- skip this step"
else
   echo "Now run FGS for variant reads..."
   $fgs -genome=$exp_d.mismatched.higherthan1.fasta -out=$exp_d.mismatched.fgs -complete=0 -train=complete -thread=$thread
   $fgs -genome=$exp_d.uncordant.fasta -out=$exp_d.uncordant.fgs -complete=0 -train=complete -thread=$thread
fi

if [ -f $exp_d.var2pep.fasta ]; then
   echo "Extraction of reads similar to identitied proteins already done -- skip this step"
else
   if [ -s $exp_d.uncordant.rapsearch.m8 ]; then
      echo "Rapsearch against identified proteins already done -- skip this step"
   else
      echo "Now run rapsearch..."
      $rap_d/prerapsearch -d $exp_d.graph2pro.fasta -n $exp_d.graph2pro.rapsearch
      $rap -q $exp_d.mismatched.fgs.faa -d $exp_d.graph2pro.rapsearch -o $exp_d.mismatched.rapsearch -z $thread > $exp_d.mismatched.errormessage
      $rap -q $exp_d.uncordant.fgs.faa -d $exp_d.graph2pro.rapsearch -o $exp_d.uncordant.rapsearch -z $thread > $exp_d.uncordant.errormessage
   fi
   echo "Now extract variant peptides..."
   python $pgm_d/pyscript/parseMismatch.py $exp_d.mismatched.rapsearch.m8 $exp_d.mismatched.fgs.faa $exp_d.mismatched.rapsearch.fasta $ide
   python $pgm_d/pyscript/parseMismatch.py $exp_d.uncordant.rapsearch.m8 $exp_d.uncordant.fgs.faa $exp_d.uncordant.rapsearch.fasta $ide


   #make var2pep database
   cat $exp_d.mismatched.rapsearch.fasta $exp_d.uncordant.rapsearch.fasta > $exp_d.var2pep.fasta
fi

if [ -s $exp_d.var2pep.mzid ]; then
   echo "MSGF+ against var2pep database already done -- skip this step"
else
   rm -f $exp_d.var2pep.mzid
   rm -f $exp_d.var2pep.revCat*
   echo "Now run MSGF+ against variants..."
   echo java -Xmx$ram -jar $msgf -s $mgf -o $exp_d.var2pep.mzid -d $exp_d.var2pep.fasta -inst 1 -t 15ppm -ti -1,2 -mod $pgm_d/MSGF+/Mods_normal.txt -ntt 1 -tda 1 -maxCharge 7 -minCharge 1 -addFeatures 1 -n 1 -thread $thread
   java -Xmx$ram -jar $msgf -s $mgf -o $exp_d.var2pep.mzid -d $exp_d.var2pep.fasta -inst 1 -t 15ppm -ti -1,2 -mod $pgm_d/MSGF+/Mods_normal.txt -ntt 1 -tda 1 -maxCharge 7 -minCharge 1 -addFeatures 1 -n 1 -thread $thread
   java -Xmx$ram -cp $msgf edu.ucsd.msjava.ui.MzIDToTsv -i $exp_d.var2pep.mzid -showDecoy 1
fi

python $pgm_d/pyscript/parseFDR_o.py $exp_d.var2pep.tsv $fdr
python $pgm_d/pyscript/parseFDR_o_peptide.py $exp_d.var2pep.tsv $fdr

#combine graph2pro and var2pep search results

python $pgm_d/pyscript/getUniquePeptides_files.py -o ${exp_d}.speFDR.final-peptide.txt -i $exp_d.graph2pro.tsv.$fdr.tsv -v $exp_d.var2pep.tsv.$fdr.tsv
python $pgm_d/pyscript/getUniquePeptides_files.py -o ${exp_d}.pepFDR.final-peptide.txt -i $exp_d.graph2pro.tsv.peptide.$fdr.tsv -v $exp_d.var2pep.tsv.peptide.$fdr.tsv
#clean up big intermediate files
cleanup="yes"
if [ "${cleanup}" = "yes" ]; then
	rm -f $exp_d.mismatched.rapsearch.m8 $exp_d.uncordant.rapsearch.m8
	rm -f $exp_d.mismatched.rapsearch.aln $exp_d.uncordant.rapsearch.aln
	rm -f $exp_d.step1.fasta.1 $exp_d.step1.fasta.2 $exp_d.step1.fasta.12
	rm -f $exp_d.sam
fi
 
echo "Unique peptides saved to ${exp_d}.final-peptide.txt"
echo "Graph2Pro-Var pipeline completed!"
