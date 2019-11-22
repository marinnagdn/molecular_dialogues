#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re
import pickle

#******************************************

def make_dict_metabolites(xml_file, dict_bacteria_kegg,dict_kegg_name):

    bacteria_name = re.split(".xml",re.split("/",xml_file)[-1])[0]

	nameHandle = open(xml_file,"r")

    list_kegg = []

	for line in nameHandle :

		if re.search("<species.*compartment=\"C_e\"", line):
			split_line = re.split('"',line)
			name = split_line[3]
        if re.search("KEGG Compound", line):
            split_line = re.split(" ",line)[2]
            kegg_name = re.split(";",split_line)[0]
        dict_kegg_name[kegg_name]=name #on écrase en partant du principe que c'est moins long d'affecter que de chercher si la valeur existe déjà, et de toute manière un id = un seul name
        list_kegg.append(kegg_name)

    dict_bacteria_kegg[bacteria_name] = list_kegg

	nameHandle.close()

	return dict_bacteria_kegg,dict_kegg_name

#********************************************

bacteria_name = re.split(".xml",re.split("/",sys.argv[1])[-1])[0]

dict_bacteria_kegg = {}
dict_kegg_name = {}

for xml in sys.argv[1]:
    make_dict_metabolites(xml,dict_bacteria_kegg,dict_kegg_name)


with open('dict_Bacteria_KEGG_exc', 'wb') as f:
	 mon_pickler = pickle.Pickler(f)
	 mon_pickler.dump(dict_bacteria_kegg)

with open('dict_KEGG_name_exc', 'wb') as f:
	 mon_pickler = pickle.Pickler(f)
	 mon_pickler.dump(dict_kegg_name)

    


