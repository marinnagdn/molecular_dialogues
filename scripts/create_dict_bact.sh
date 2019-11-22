
#!/usr/bin/env bash

# No args needed
# EXAMPLE OF COMMAND LINE:
# From your results file with all interactions list make :
# ./get_molecular_dialogues.sh

IFS=$'\n'  #to read whole line and not cut line by white spaces
for i in $(cat list_models.txt); do
    path=$(find Database/Embl -name "$i")
	gunzip -c $path > ./$(basename $i .gz) 2> /dev/null
done
cut -d"." -f 1,2 list_models.txt > list_models_xml.txt
python3 create_dic_bact_KEGG.py list_models_xml.txt
