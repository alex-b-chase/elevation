#!/bin/bash

BASE=/Volumes/Bennett_BACKUP/Research/curtobacterium/curto-elevation
OUTDIR=$BASE/patric_genomes
REFDB=/Volumes/Bennett_BACKUP/Research/reference_db
ESSEN=$REFDB/essen_singlecopy

mkdir -p $OUTDIR
cd $OUTDIR

rm -rf $OUTDIR/prodigal_annotation
rm -rf $OUTDIR/markergenes
# rm -rf $OUTDIR/genomes
rm -rf $OUTDIR/ftp.patricbrc.org
mkdir -p $OUTDIR/genomes
mkdir -p $OUTDIR/prodigal_annotation
mkdir -p $OUTDIR/markergenes

# add "rep genome" designation to each line to match format of other files
# only grab columns that we want and reorder them to match the other files
# remove lines that are blank in certain fields | remove double quotes | change column headers
# finally, take out genomes from common genera that we already have processed!
# remove lines with identical genome names in column 2
cat $BASE/PATRIC_genome_rep.txt | awk -v OFS='\t' '{print $0, "1"}' | \
awk -v OFS='\t' -F'\t' '{print $1,$2,$69,$4,$5,$16,$17,$18,$19,$29,$30,$31,$32,$36,$39,$65}' | awk -F'\t' '$2' | \
sed 's/\"//g' | sed -e '1s/1/rep_genome/' | grep -v "Curtobacterium" | \
tr -d '-' | awk -F"\t" '!seen[$2]++' > $BASE/patric_rep_genomes.txt

# FINALLY, have a unique set of genomes to work with

# cut -f1 $BASE/patric_rep_genomes.txt | tail +2 | sort | uniq > download.txt


# # download the files from PATRIC
# while read line
# do
# 	filename=$line.fna

# 	cd $OUTDIR/genomes

# 	# if we already downloaded the file, skip it
# 	if [ ! -f "$filename" ]
# 	then
# 		cd $OUTDIR
# 		wget -cqNr -t 45 -A "*.fna" "ftp://ftp.patricbrc.org/patric2/genomes/$line/" -P .
# 		echo "done with $line"

# 	else
# 		echo "done with $line"
# 		continue
# 	fi

# done < download.txt

# rm -f download.txt
# cd $OUTDIR/ftp.patricbrc.org

# # move downloaded genomes to new directory
# find . -type f -name '*.fna' -exec mv {} $OUTDIR/genomes \;

# cd $OUTDIR/genomes
# rm -r $OUTDIR/ftp.patricbrc.org

#################################################################################
####################### 	CURTOBACTERIUM GENOMES 	#############################
#################################################################################

# these have been previously analyzed thoroughly and picked one rep genome
# also, included Kristin's isolates that Sydney sequenced

genomeMAP=$BASE/mapping.genomeFile.txt

CURTODIR=$BASE/survey_isolates/prodigal_annotation
LOMADIR=$BASE/sydney_isolates


# subset only the representative genomes from Curtobacterium
awk -F $'\t' 'BEGIN {OFS = FS} {if ($5 == "1") print $0;}' $genomeMAP > $BASE/curtoLR_rep_genomes.txt


while read newname filename genus directory rep
do

	if [[ "${directory}" == "curto" ]] 
		then
		cd $CURTODIR 
		cp $filename $OUTDIR/genomes/$newname.fna
	else 
		cd $LOMADIR
		cp $filename $OUTDIR/genomes/$newname.fna
	fi

done < $BASE/curtoLR_rep_genomes.txt

rm $BASE/curtoLR_rep_genomes.txt


#################################################################################
#################################################################################

cd $OUTDIR/genomes

# annotate and screen for marker genes
refmgs=$ESSEN/Essen_SingleCopy.txt
cat $refmgs | awk '{if ($6 == "1") print $0;}' > $OUTDIR/prodigal_annotation/markers.sub.txt

mgs=$OUTDIR/prodigal_annotation/markers.sub.txt
check=$(cat $mgs | wc -l)


# subset for naming purposes only
cut -f1-2 $BASE/patric_rep_genomes.txt | tail +2 > $OUTDIR/patric_list.txt
genomeMAP=$OUTDIR/patric_list.txt


genomecount=$(find . -type f -name '*.fna'  | wc -l)
echo "DONE! There are ${genomecount} genomes ready to be annotated!"


# remove old files
rm -f *temp*
rm -f $OUTDIR/deleted_genomes.txt
rm -f $OUTDIR/excess_markers.txt
rm -f $OUTDIR/duplicate_markers.txt
rm -f $OUTDIR/final_passed_genomes.txt
rm -f *.markers.faa
rm -f *ids.txt
rm -f *.faa

# start the annotation and screening!

for f in *.fna
do
	initial_start=$(date +%s)

	base=${f%.fna}
	patric=$(echo $f | cut -f1 -d'.')
	# check if genome is PATRIC or not - PATRIC genomes start with numbers
	if ! [[ $patric =~ ^[0-9]+$ ]]
	then
		output=$base
	else
		output=$(cat $genomeMAP | grep -w "$base" | cut -f2 | tr -d "[-/%,;\(\):=\.\\\[]\"\']" | sed 's/  */_/g')

		# if we can't find the $output variable in the mappingFile, output that for reference - should be nothing!
		if [ -z "$output" ] 
		then
			echo $base
			echo $output
			echo ''
		else
			:
		fi
	fi

	# # translate using prodigal, skip empty files
	if [[ ! -s $f ]]
		then
			echo $f >> $OUTDIR/somethingwrong.txt
			continue
		else
			echo "Annotating ${output}..."

			prodigal -i $f -a $output.faa -f gff -q > $output.gff

	fi

	# search against the Banfield marker gene database
	hmmsearch --tblout $output.temp.txt -E 1e-10 --cpu 4 $ESSEN/Essen_SingleCopy.hmm $output.faa > /dev/null

	# now we can subset by the phylomarkers
	while read acc bla blah blah blahh blahhh protein
	do
		cat $output.temp.txt | grep $acc | cut -f1 -d' ' | sort | uniq >> $output.temp2.txt
		cat $output.temp.txt | grep $acc | tr -s ' ' | cut -f3 -d' ' | sort | uniq >> $output.temp3.txt

		filename=$(echo $output.faa)
		sequence=$(cat $output.temp.txt | grep $acc | cut -f1 -d' ' | sort | uniq)
		prot=$(echo $protein | cut -f2-3 -d' ' | sed 's/ /_/g' | cut -f1 -d'-')
		echo -e $filename'\t'$sequence >> ${prot}p.ids.txt

	done < $mgs

	# check to make sure genomes only have one copy of each...
	markercheck=$(cat $output.temp2.txt | wc -l)
	protcheck=$(cat $output.temp3.txt | wc -l)

	if [ $markercheck -lt $check ]
	then
		echo -e $f'\t'$output'\t'$markercheck >> $OUTDIR/deleted_genomes.txt
		rm *${output}*
		continue

	elif [ $markercheck -gt $check ]
	then
		echo -e $f'\t'$output'\t'$markercheck >> $OUTDIR/excess_markers.txt
		rm *${output}*
		continue

	elif [ $markercheck -ne $protcheck ]
	then
		echo -e $f'\t'$output'\t'$protcheck >> $OUTDIR/duplicate_markers.txt
		rm *${output}*
		continue		

	else
		echo "Genome ${output} passed!"
		cat $BASE/patric_rep_genomes.txt | grep -w "$base" >> $OUTDIR/final_passed_genomes.txt
		

	fi		

	# filter out the marker genes from the genome
	search-fasta.py -i $output.faa -m $output.temp2.txt -o $output.markers.faa


	total_end=$(date +%s)
	total_runtime=$(echo "$total_end - $initial_start" | bc -l)

	echo "Done annotating $output - total time $total_runtime seconds"

	rm *.temp*
	rm *.gff

done

# subset the marker gene annotations to separate folder
rm -rf genome_markers
mkdir -p genome_markers
find . -name "*.markers.faa" -maxdepth 1 -exec mv {} genome_markers/ \;

# move the genomes and marker gene annotations to own file
find . -name "*.faa" -maxdepth 1 -exec mv {} $OUTDIR/prodigal_annotation \;
find . -name "*.ids.txt" -maxdepth 1 -exec cp {} $OUTDIR/prodigal_annotation \;


################################################################################################
# go through the reference ids.txt files to find the genome file and sequence corresponding
# to the AA sequence of the marker genes
# output the results into a marker gene specific AA file with all seqs in there
################################################################################################

cd $OUTDIR/prodigal_annotation

genomecount=$(find . -type f -name '*.faa'  | wc -l)
echo "DONE! There are ${genomecount} genomes that were downloaded and processed from PATRIC"

# remove all old files that may have previously been ran
rm -f *.ids.faa

# loop through the lines with the hit marker genes and subset by the marker genes
for f in *.ids.txt
do
	filename=${f%.ids.txt}

	echo "Processing marker gene $filename"

	while read gen sequence 
	do
		genome=${gen%.faa}
		sequence2=$(echo $sequence | cut -f1 -d' ')
		#echo $sequence2

		# some genomes were deleted 
		# so if genome was deleted, do nothing (:)
		if [ -f $gen ]
			then 
				awk -v "seq=$sequence2" -v RS='>' '$1 == seq {print RS $0}' $gen | \
				sed "s/[>].*$/>$genome/" >> $filename.ids.faa
			else
				:
		fi
	done < $f
done

# remove all old folders and files related to this
# move the newly produced AA files to a new folder

mv *.ids.faa $OUTDIR/markergenes
cp *.ids.txt $OUTDIR/markergenes

################################################################################################
# align all the individual protein files for tree construction
################################################################################################


cd $OUTDIR/markergenes

rm -f *.aln
rm -f *-no-dups.fasta


for f in *.ids.faa
do
	newname="${f/%.ids.faa/.total.fa}"
	alignment="${f/%.ids.faa/.total.aln}"

	# remove all weird characters (RAXML does not like them!)
	cat $f | tr -d "[ -%,;\(\):=\.\\\[]\"\']" | sed "s/\*//g" > "$newname"

	# align the sequences to create aligned fasta files
	clustalo -i $newname -o $alignment

	# remove duplicate IDs just in case there were any relic mistakes
	no_dups.sh $alignment
	# this will output our final file == "*-no-dups.fasta"

	rm -f $f
done



################################################################################################
# build a consensus reference tree by combining all 30 marker genes into one file
# construct the tree
################################################################################################

# combine the fasta files together for tree construction keeping the proteins correctly ordered
# all individual alignment files must have the same header

catfasta2phyml.pl -f *.fasta  > concat.aligned.filtered.fa

# build the tree using the aligned file and output 100 bootstrap replications
# kept crashing, running on the other lab computer - WORKED so run separately!

rm -f RAxML*
rm -rf $OUTDIR/markergenes/master_tree/
mkdir $OUTDIR/markergenes/master_tree/

mv concat.aligned.filtered.fa $OUTDIR/markergenes/master_tree/
cd $OUTDIR/markergenes/master_tree/

FastTree concat.aligned.filtered.fa > concat.aligned.filtered.tre 

# try and run on HPC, no way will this work on desktop
# raxmlHPC-PTHREADS -s concat.aligned.filtered.fa -m PROTGAMMABLOSUM62 -n markergene -x 100 -# 100 -p 4321 -f a -T 6








