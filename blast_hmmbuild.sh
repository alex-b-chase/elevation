#!/bin/bash

BASE=/Volumes/Bennett_BACKUP/Research/curtobacterium/curto-elevation/patric_genomes/
REFDIR=$BASE/markergenes
OUTDIR=$BASE/referenceDB 

mkdir -p $OUTDIR

# use the final concat file to subset genomes and grab taxIDs from mappingFile
cd $REFDIR/master_tree

print-fasta-id.py concat.aligned.filtered.fa

# double sure it matches up
cat concat.aligned.filtered.fa | grep '>' | wc -l
echo "should be 5433 genomes"
echo ""

refgen=$BASE/ID2names.fix.txt
outgen=$OUTDIR/final.genomes.txt

# rm -f $outgen

# # need to do exact string match to avoid similar ones (i.e. CAG57 grep both CAG57 and CAG579)
# while read line
# do
# 	count=`cat $refgen | grep "\b${line}\b" | wc -l`
	
# 	if [ $count -eq 1 ]
# 	then
#     	cat $refgen | grep -w $line >> $outgen
#     elif [ $count -gt 1 ]
#     then
#     	echo "${line} is repetitive. Check this out yourself..."
# 	else
#     	#echo "${line} not found. Look up taxID on NCBI and add manually. Sort ${outgen} by first column and replace 0 with taxID"
#     	echo "0:${line}" | tr ':' '\t' >> $outgen
# 	fi

# done < concat.aligned.filtered_ids.txt

# now, if you are done with that and are satisfied, can move on and comment above out...
# save the resulting final file as $OUTDIR/allgenomeIDs.txt
# we will use this to make the BLAST DB mapping file and rename sequences to make into BLAST DBs

wc -l $OUTDIR/allgenomeIDs.txt
refgen=$OUTDIR/allgenomeIDs.txt

# make newnames for each fasta file with new headers

cd $REFDIR
rm -f $REFDIR/newnames.txt
rm -f $OUTDIR/FINAL_markers.map

# format for newnames == "gnl|BACT|$taxID[i]_$protein taxon=$taxID, $genomeID"


for f in *.good.final.fa
do
	print-fasta-id.py $f 

	protein=${f%.good.final.fa}

	count=1
	rm -f temp.txt

	while read line
	do
		refgenome=$(cat $refgen | grep "\b${line}\b")

		genomeID=$(echo ${refgenome} | cut -f2 -d' ')
		taxID=$(echo ${refgenome} | cut -f1 -d' ')

		echo "gnl|BACT|${taxID}.${count}_${protein}:taxon=${taxID},:${genomeID}" >> temp.txt
		echo "gnl|BACT|${taxID}.${count}_${protein}:${taxID}" | tr ':' '\t' >> $OUTDIR/FINAL_markers.map

		count=`expr $count + 1`

	done < $protein.good.final_ids.txt

	paste $protein.good.final_ids.txt temp.txt > newnames.$protein.txt

	# now rename the fasta files and output for BLAST DB
	fasta-rename.py $f newnames.$protein.txt $protein.blast.temp.fasta 
	cat $protein.blast.temp.fasta | tr ':' ' ' > $OUTDIR/$protein.blast.fasta 

	# rm -f temp.txt
	# rm -f newnames.$protein.txt
	# rm -f $protein.blast.temp.fasta

done


cd $OUTDIR

cat *.blast.fasta > total_markers.faa
rm *.blast.fasta 

# create a custom local database
makeblastdb -in total_markers.faa -dbtype 'prot' -out total_markergene -parse_seqids -taxid_map FINAL_markers.map

# check to make sure it worked
blastdbcmd -db total_markergene -entry all -outfmt "%T"

# should get a long list of tax ID


# now build HMMer profiles for DB as well...

cd $REFDIR/master_tree

for f in *.good.final.aln 
do
	cd $REFDIR/master_tree
	protein=${f%.good.final.aln}

	cp $f $OUTDIR
	cd $OUTDIR

	hmmbuild $protein.hmm $f >> $protein.log

	rm $f 

done


