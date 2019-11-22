#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re
import pickle

#******************************************

def make_dict_metabolites(xml_file):

	dict_metabolites = {}

	nameHandle = open(xml_file,"r")

	for line in nameHandle :

		if re.search("<species.*name=", line):
			split_line = re.split('"',line)
			id = split_line[1]
			name = split_line[9]
			dict_metabolites [id] = name

	nameHandle.close()

	return dict_metabolites

#********************************************

def make_dict_reactions(xml_file):

	dict_reactions = {}

	nameHandle = open(xml_file,"r")

	for line in nameHandle :

		if re.search("<reaction.*id=", line):
			split_line = re.split('"',line)
			id = split_line[1]
			name = split_line[7]
			dict_reactions[id] = [name]

	nameHandle.close()

	return dict_reactions

#********************************************

metabolites_human = make_dict_metabolites(sys.argv[1])

reactions_human = make_dict_reactions(sys.argv[1])

with open('Homo_Sapiens_metabolites_database', 'wb') as f:
	 mon_pickler = pickle.Pickler(f)
	 mon_pickler.dump(metabolites_human)

with open('Homo_Sapiens_reactions_database', 'wb') as f:
	 mon_pickler = pickle.Pickler(f)
	 mon_pickler.dump(reactions_human)