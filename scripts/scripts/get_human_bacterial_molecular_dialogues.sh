#!/usr/bin/env bash

# ARGUMENTS :
# Give the tsv files of interest as arguments
# EXEMPLE OF COMMAND LINE :
# ./get_human_bacterial_molecular_dialogues.sh ../tsv_files/NW.tsv ../tsv_files/ANT1.tsv ../tsv_files/ANT2.tsv

IFS=$'\n'
cpt=0

mkdir RESULTS_get_human_bacterial_molecular_dialogues 2> /dev/null

for tsv in "$@"; do
	nameHandle=$(basename $tsv .tsv)
	awk 'FNR==1{print $1;next} {for(i=2;i<=NF;i++){sum+=$i};print $1,sum;sum=""}'  $tsv > ${nameHandle}_list_bacterias.txt #sum of abundances by bacteria
	# awk '$2 > 0' ${nameHandle}_list_bacterias.txt > tmp && mv tmp ${nameHandle}_list_bacterias.txt #get only at least once present bacterias in one individual within the sample
		awk -F "g__" '{print $2}'  ${nameHandle}_list_bacterias.txt | awk -F ";" '{print $1}' > ${nameHandle}_genre.txt #extract genre
		awk -F "s__" '{print $2}'  ${nameHandle}_list_bacterias.txt | awk -F ";" '{print $1}' > ${nameHandle}_species.txt #extract species
		paste ${nameHandle}_genre.txt ${nameHandle}_species.txt > ${nameHandle}_genre_species.txt #merge genre and species
		paste ${nameHandle}_list_bacterias.txt ${nameHandle}_genre_species.txt | awk {'print $3 " " $4 "\t" $2'} > tmp && mv tmp ${nameHandle}_list_bacterias.txt #make new list files of bacterias with Genre Species ID at the first column
		rm ${nameHandle}_genre_species.txt ${nameHandle}_genre.txt ${nameHandle}_species.txt #clean up
		cpt=$(($cpt + 2))
done

#once we have our lists of bacterias by subgroup 

#create results list file
for tsv in "$@"; do
	nameHandle=$(basename $tsv .tsv)
	output="${output}_${nameHandle}"
	subgroups="${subgroups} $nameHandle"
done
output="$(echo $output | cut -d _ -f 2-).txt"

echo -en "Bacteria\t" > $output
for tsv in "$@"; do
	nameHandle=$(basename $tsv .tsv)
	printf "$nameHandle\t" >> $output
done

paste *_list_bacterias.txt* > tmp_list_bacterias.txt

#JUSTE POUR NOS DONNEES !!!
#car j'arrive pas à généraliser
awk  '{print $1 " " $2 "\t"  $3 "\t" $6 "\t" $9}' tmp_list_bacterias.txt >> $output

rm tmp_list_bacterias.txt

tr -d "[]" < $output > tmp && mv tmp $output #remove crochets

for bacteria in $(awk ' NR > 1 {print $1"_"$2}' $output);do
	echo $bacteria
	path=$(find ../Database/Embl -name "*$bacteria*" | head -n 1 -q) #only first match, pas les autres souches #avec un _ entre le genre et l'espèce
	if [  -z $path1 ]; then #found in the database

		bacteria_prefix_file=$(basename $path .gz) # prefix file

		gunzip -c $path > ${bacteria_prefix_file} 2> /dev/null # decompress

		python3 ../scripts/bin/get_human_bacterial_molecular_dialogues.py  ${bacteria_prefix_file} ../Database/Recon3D/Homo_sapiens.xml RESULTS_get_human_bacterial_molecular_dialogues/

		rm ${bacteria_prefix_file} 2> /dev/null 
	
	else 

		echo "$bacteria not found in the database"

	fi

done

