#!/usr/bin/env bash

# No args needed
# EXAMPLE OF COMMAND LINE:
# From your results file with all interactions list make :
# ./get_molecular_dialogues.sh

IFS=$'\n'  #to read whole line and not cut line by white spaces


# for i in unique_interactions.txt; do 
for i in $(ls ../Results/*interactions.txt); do 
	#NR > 1 : skip first line 
	echo -e "----------- $i : -----------"
	dossier=$(basename $i .txt)

	mkdir $dossier  2> /dev/null
	# awk 'NR>1 {print $1,$2}' $i > $dossier/left; #left bacterias
	# awk 'NR>1 {print $3,$4}' $i > $dossier/right; #right bacterias

	for line in $(cat $i | tail -n +2); do #tail -n +2 to ignore header first line
		#on récupère les noms des bactéries dans des variables
		# echo $line

		bacteria1=$(echo $line | cut -f 1); bacteria2=$(echo $line | cut -f 2)
		echo "$bacteria1 - $bacteria2"
		genre1=$(echo $bacteria1 | cut -f1 -d " "); species1=$(echo $bacteria1 | cut -f2 -d " ")
		genre2=$(echo $bacteria2 | cut -f1 -d " "); species2=$(echo $bacteria2 | cut -f2 -d " ")

		# on récupère les chemins des modèle au format xml.gz
		path1=$(find ../Database/Embl -name "$genre1\_$species1*" | head -n 1 -q) #only first match, pas les autres souches #avec un _ entre le genre et l'espèce

		path2=$(find ../Database/Embl -name "$genre2\_$species2*" | head -n 1 -q) #only first match, pas les autres souches

		if [ ! -z $path1 ] && [ ! -z $path2 ]; then

			bacteria1_prefix_file=$(basename $path1 .gz) #decompress prefix file
			bacteria2_prefix_file=$(basename $path2 .gz)

			#decompress
			gunzip -c $path1 > ${dossier}/${bacteria1_prefix_file} 2> /dev/null 
			gunzip -c $path2 > ${dossier}/${bacteria2_prefix_file} 2> /dev/null

			cd $dossier

			# echo ${bacteria1_prefix_file}
			# echo ${bacteria2_prefix_file}

			python3 bin/get_bacterial_molecular_dialogues.py ${bacteria1_prefix_file} ${bacteria2_prefix_file} $(pwd)

			# for file in $(ls *.list); do
			# 	for m in $(cat $file); do
			# 		grep "id=\"${m}\"" $bacteria1_prefix_file | cut -d \" -f 4  >> $(basename $file .list)_metabolites.txt #obtenir le nom complet du métabolite par rapport à l'identifiant
			# 		sort $(basename $file .list)_metabolites.txt | uniq >> tmp && mv tmp $(basename $file .list)_metabolites.txt 
			# 		# paste $file $(basename $file .list)_metabolites.txt
			# 	done
			# done

			#clean up
			rm ${bacteria1_prefix_file} 2> /dev/null 
			rm ${bacteria2_prefix_file} 2> /dev/null
			# rm *.list 2> /dev/null
		
			cd .. #come back to list_interactions folder

		else 
			if [[ -z "$path1" ]];then
				echo "* $bacteria1 file not found in the database"
			fi
			if [[ -z "$path2" ]];then
				echo "* $bacteria2 file not found in the database"
			fi
		fi

	done

done

