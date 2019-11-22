#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#to get cytosol -> extracellular and extracellular -> cytosol metabolites from a xml metabolism bacteria file

# python3 extract_c_e_metabolites.py /home/marinna/M2/Bioinfo_appliquée_2/embl_gems-master/models/a/acidobacterium/toy_Acidobacterium_ailaaui_PMMR2.xml

import sys
import re
import pickle

#*************************************

def make_dict_metabolites(xml_file, dict_human_metabolites):

	nameHandle = open(xml_file,"r")

	dict_metabolites = dict_human_metabolites

	for line in nameHandle :

		if re.search("<species.*name=", line):
			split_line = re.split('"',line)
			id = split_line[1]
			if id not in dict_metabolites.keys():
				name = split_line[3]
				dict_metabolites [id] = name

	nameHandle.close()

	return dict_metabolites

#*************************************

def make_dict_reactions(xml_file, dict_human_reactions):

	dict_reactions = dict_human_reactions

	nameHandle = open(xml_file,"r")

	for line in nameHandle :

		if re.search("<reaction.*id=", line):
			split_line = re.split('"',line)
			id = split_line[1]
			name = split_line[3]
			if id not in dict_reactions.keys():
				dict_reactions[id] = [name]
			else: 
				dict_reactions[id].extend([name])

	nameHandle.close()

	return dict_reactions


#*************************************

def get_dict_reactions(xml_file, metabo_bacteria):

	has_found_reaction = False
	has_found_reactant = False
	has_found_product = False

	is_reversible = False

	description_reaction = []
	description_reactants = []
	description_products = []

	#reversible
	list_extracellular_M = []
	# list_periplasm_M = []
	# list_cytosol_M = []

	#non reversible
	list_reactants_extracellular_M = []
	list_reactants_periplasm_M = []
	list_reactants_cytosol_M = []
	list_products_extracellular_M =[]
	list_products_periplasm_M = []
	list_products_cytosol_M = []

	dict_reactions = {} #dict with id reactions as keys, and list of lists 1) reactants 2) products as values
	dict_detailed_reactions = {}

	nameHandle = open(xml_file,"r")

	for line in nameHandle :
		if re.search(r"</listOfReactionsQ>",line): #end of reactions
			break
		elif re.search(r"<reaction id=",line): #new reaction
			reaction_id = re.split('"',line)[1]	
			reaction_name = re.split('"',line)[3]	

			#Reversible reaction ???
			if re.search("reversible=\"true\"",line):
					is_reversible = True

			if not re.search("EX",reaction_id): #delete EX reactions
				has_found_reaction = True #lets get informations of this reaction
		elif has_found_reaction:
			description_reaction.append(line) #accumule informations of the reactions

		#End of the current reaction : let's analyse it
		if re.search(r"</reaction>",line): 

			if is_reversible: #pourquoi j'ai utilisé des dictionnaires déjà là ??????????????????

				for line2 in description_reaction:

					m = re.findall(r"M_.*_e\"",line2)
					if m != []: #if found metabolite
						list_extracellular_M.extend(m[:-1])

					# m = re.findall(r"M_.*_c\"",line2)
					# if m != []: #if found metabolite
					# 	list_cytosol_M.extend(m)

					# m = re.findall(r"M_.*_p\"",line2)
					# if m != []: #if found metabolite
					# 	list_periplasm_M.extend(m)

				# Add in dictionaries for reversible : only compartments known		

				# dict_reactions[reaction_id] = [list(dict_extracellular_M.keys()),list(dict_cytosol_M.keys()), list(dict_periplasm_M.keys())]

				# dict_detailed_reactions[reaction_name] = [[metabo_bacteria[x] for x in list(dict_extracellular_M.keys())],  [metabo_bacteria[x] for x in list(dict_cytosol_M.keys())], [metabo_bacteria[x] for x in list(dict_periplasm_M.keys())]]

				dict_reactions[reaction_id] = list_extracellular_M

				dict_detailed_reactions[reaction_name] = [metabo_bacteria[x] for x in list_extracellular_M]

				description_reaction = [] #erase description
				has_found_reaction = False
				is_reversible = False

				list_extracellular_M = []
				list_periplasm_M = []
				list_cytosol_M = []

			else:

				for line3 in description_reaction:

					#Reactants of the reaction
					if re.search("<listOfReactants>",line3):
						has_found_reactant = True
					elif has_found_reactant:
						description_reactants.append(line3)

					if re.search(r"</listOfReactants>",line3):
						
						for reactant in description_reactants:

							r = re.findall(r"M_.*_p\"",reactant)
							if r != []:
								for r2 in r:
									list_reactants_periplasm_M.append(r2[:-1])

							r = re.findall(r"M_.*_c\"",reactant)
							if r != []:
								for r2 in r:
									list_reactants_cytosol_M.append(r2[:-1])


							r = re.findall(r"M_.*_e\"",reactant)
							if r != []:
								for r2 in r:
									list_reactants_extracellular_M.append(r2[:-1])


						has_found_reactant = False
						description_reactants = []

					#Products of the reaction
					if re.search(r"<listOfProducts>",line3):
						has_found_product = True
					elif has_found_product:
						description_products.append(line3)

					if re.search(r"</listOfProducts>",line3):

						for product in description_products:

							p = re.findall(r"M_.*_p\"",product)
							if p != []:
								for p2 in p:
									list_products_periplasm_M.append(p2[:-1])

							p = re.findall(r"M_.*_c\"",product)
							if p != []:
								for p2 in p:
									list_products_cytosol_M.append(p2[:-1])

							p = re.findall(r"M_.*_e\"",product)
							if p != []:
								for p2 in p:
									list_products_extracellular_M.append(p2[:-1])

						has_found_product = False
						description_products = []

				############################

				# Analysis all metabolite lists.

				# ID to NAME Metabolites
				reactants_extracellular_M = [metabo_bacteria[x] for x in list_reactants_extracellular_M]
				reactants_periplasm_M = [metabo_bacteria[x] for x in list_reactants_periplasm_M]
				reactants_cytosol_M = [metabo_bacteria[x] for x in list_reactants_cytosol_M]
				products_extracellular_M = [metabo_bacteria[x] for x in list_products_extracellular_M]
				products_periplasm_M = [metabo_bacteria[x] for x in list_products_periplasm_M]
				products_cytosol_M = [metabo_bacteria[x] for x in list_products_cytosol_M]

				#ENTERING metabolites
				if (len(list_reactants_extracellular_M) != 0):
					# in_M.extend(list_reactants_extracellular_M)
					# print(reaction_id)
					dict_reactions[reaction_id] = [list_reactants_extracellular_M,list_reactants_periplasm_M,list_reactants_cytosol_M],[list_products_extracellular_M,	list_products_periplasm_M,	list_products_cytosol_M]

					dict_detailed_reactions[reaction_name] = [reactants_extracellular_M,reactants_periplasm_M,reactants_cytosol_M],[products_extracellular_M,	products_periplasm_M,	products_cytosol_M]

				#EXITING metabolites
				if (len(list_products_extracellular_M)  != 0):

					dict_reactions[reaction_id] = [list_reactants_extracellular_M,list_reactants_periplasm_M,list_reactants_cytosol_M],[list_products_extracellular_M,	list_products_periplasm_M,	list_products_cytosol_M]

					dict_detailed_reactions[reaction_name] = [reactants_extracellular_M,reactants_periplasm_M,reactants_cytosol_M],[products_extracellular_M,	products_periplasm_M,	products_cytosol_M]

				#Empty lists
				list_reactants_extracellular_M = []
				list_reactants_periplasm_M = []
				list_reactants_cytosol_M = []
				list_products_extracellular_M = []
				list_products_periplasm_M = []
				list_products_cytosol_M = []
				reactants_extracellular_M = []
				reactants_periplasm_M = []
				reactants_cytosol_M = []
				products_extracellular_M = []
				products_periplasm_M = []
				products_cytosol_M = []

				description_reaction = [] #erase description
				has_found_reaction = False

	nameHandle.close()

	return dict_reactions, dict_detailed_reactions

#*************************************

def get_dict_reactions_homo_sapiens(xml_file, metabo_bacteria,reactions_bacteria):

	has_found_reaction = False
	has_found_reactant = False
	has_found_product = False

	is_reversible = False

	description_reaction = []
	description_reactants = []
	description_products = []

	#reversible
	list_extracellular_M = []
	# list_periplasm_M = []
	# list_cytosol_M = []

	#non reversible
	list_reactants_extracellular_M = []
	list_reactants_periplasm_M = []
	list_reactants_cytosol_M = []
	list_products_extracellular_M =[]
	list_products_periplasm_M = []
	list_products_cytosol_M = []

	dict_reactions = {} #dict with id reactions as keys, and list of lists 1) reactants 2) products as values
	dict_detailed_reactions = {}

	nameHandle = open(xml_file,"r")

	for line in nameHandle :
		if re.search(r"</listOfReactionsQ>",line): #end of reactions
			break
		elif re.search(r"<reaction id=",line): #new reaction
			reaction_id = re.split('"',line)[1]	
			reaction_name = re.split('"',line)[7]

			#Reversible reaction ???
			if re.search("reversible=\"true\"",line):
					is_reversible = True	

			if not re.search("EX",reaction_id): #delete EX reactions
				has_found_reaction = True #lets get informations of this reaction
		elif has_found_reaction:
			# print (line)
			description_reaction.append(line) #accumule informations of the reactions

		#End of the current reaction : let's analyse it
		if re.search(r"</reaction>",line): 
			
			if is_reversible: 

				for line2 in description_reaction:

					m = re.findall(r"M_.*_e\"",line2)
					if m != []: #if found metabolite
						list_extracellular_M.extend(m[:-1])

					# m = re.findall(r"M_.*_c\"",line2)
					# if m != []: #if found metabolite
					# 	list_cytosol_M.extend(m)

					# m = re.findall(r"M_.*_p\"",line2)
					# if m != []: #if found metabolite
					# 	list_periplasm_M.extend(m)

				# Add in dictionaries for reversible : only compartments known		

				# dict_reactions[reaction_id] = [list(dict_extracellular_M.keys()),list(dict_cytosol_M.keys()), list(dict_periplasm_M.keys())]

				# dict_detailed_reactions[reaction_name] = [[metabo_bacteria[x] for x in list(dict_extracellular_M.keys())],  [metabo_bacteria[x] for x in list(dict_cytosol_M.keys())], [metabo_bacteria[x] for x in list(dict_periplasm_M.keys())]]

				dict_reactions[reaction_id] = list_extracellular_M

				dict_detailed_reactions[reaction_name] = [metabo_bacteria[x] for x in list_extracellular_M]

				description_reaction = [] #erase description
				has_found_reaction = False
				is_reversible = False

				list_extracellular_M = []
				list_periplasm_M = []
				list_cytosol_M = []

			else:

				if reaction_name == "Potassium transport out via proton antiport":
						print("Le potassium n'est pas réversible !!!")

				for line3 in description_reaction:

					#Reactants of the reaction
					if re.search("<listOfReactants>",line3):
						has_found_reactant = True
					elif has_found_reactant:
						description_reactants.append(line3)

					if re.search(r"</listOfReactants>",line3):
						
						for reactant in description_reactants:

							r = re.findall(r"M_.*_p\"",reactant)
							if r != []:
								for r2 in r:
									list_reactants_periplasm_M.append(r2[:-1])

							r = re.findall(r"M_.*_c\"",reactant)
							if r != []:
								for r2 in r:
									list_reactants_cytosol_M.append(r2[:-1])


							r = re.findall(r"M_.*_e\"",reactant)
							if r != []:
								for r2 in r:
									list_reactants_extracellular_M.append(r2[:-1])


						has_found_reactant = False
						description_reactants = []

					#Products of the reaction
					if re.search(r"<listOfProducts>",line3):
						has_found_product = True
					elif has_found_product:
						description_products.append(line3)

					if re.search(r"</listOfProducts>",line3):

						for product in description_products:

							p = re.findall(r"M_.*_p\"",product)
							if p != []:
								for p2 in p:
									list_products_periplasm_M.append(p2[:-1])

							p = re.findall(r"M_.*_c\"",product)
							if p != []:
								for p2 in p:
									list_products_cytosol_M.append(p2[:-1])

							p = re.findall(r"M_.*_e\"",product)
							if p != []:
								for p2 in p:
									list_products_extracellular_M.append(p2[:-1])

						has_found_product = False
						description_products = []

				############################

				# Analysis all metabolite lists.

				# ID to NAME Metabolites
				reactants_extracellular_M = [metabo_bacteria[x] for x in list_reactants_extracellular_M]
				reactants_periplasm_M = [metabo_bacteria[x] for x in list_reactants_periplasm_M]
				reactants_cytosol_M = [metabo_bacteria[x] for x in list_reactants_cytosol_M]
				products_extracellular_M = [metabo_bacteria[x] for x in list_products_extracellular_M]
				products_periplasm_M = [metabo_bacteria[x] for x in list_products_periplasm_M]
				products_cytosol_M = [metabo_bacteria[x] for x in list_products_cytosol_M]

				#ENTERING metabolites
				if (len(list_reactants_extracellular_M) != 0):
					# in_M.extend(list_reactants_extracellular_M)
					# print(reaction_id)
					dict_reactions[reaction_id] = [list_reactants_extracellular_M,list_reactants_periplasm_M,list_reactants_cytosol_M],[list_products_extracellular_M,	list_products_periplasm_M,	list_products_cytosol_M]

					dict_detailed_reactions[reaction_name] = [reactants_extracellular_M,reactants_periplasm_M,reactants_cytosol_M],[products_extracellular_M,	products_periplasm_M,	products_cytosol_M]

				#EXITING metabolites
				if (len(list_products_extracellular_M)  != 0):

					dict_reactions[reaction_id] = [list_reactants_extracellular_M,list_reactants_periplasm_M,list_reactants_cytosol_M],[list_products_extracellular_M,	list_products_periplasm_M,	list_products_cytosol_M]

					dict_detailed_reactions[reaction_name] = [reactants_extracellular_M,reactants_periplasm_M,reactants_cytosol_M],[products_extracellular_M,	products_periplasm_M,	products_cytosol_M]

				#Empty lists
				list_reactants_extracellular_M = []
				list_reactants_periplasm_M = []
				list_reactants_cytosol_M = []
				list_products_extracellular_M = []
				list_products_periplasm_M = []
				list_products_cytosol_M = []
				reactants_extracellular_M = []
				reactants_periplasm_M = []
				reactants_cytosol_M = []
				products_extracellular_M = []
				products_periplasm_M = []
				products_cytosol_M = []

				description_reaction = [] #erase description
				has_found_reaction = False

	nameHandle.close()

	return dict_reactions, dict_detailed_reactions

#***************************************

def pretty_ONE_reaction(dr, dict_detailed_reactions, id_reaction):

	one_reaction = ""

	name_reaction = dr[id_reaction][0] # [0] in case of two names for same id (list of names for a id)

	if dict_detailed_reactions[name_reaction]!=[]:

		reactants = dict_detailed_reactions[name_reaction][0]

		reactants = [x for x in reactants if x != []] #remove empty lists

		nb_reactants = len(reactants)
		cpt = 1
		for r in  reactants:
			if cpt < nb_reactants:
				one_reaction += "".join(r) + " + "
			else:
				one_reaction += "".join(r)
			cpt += 1

		one_reaction += " -> "

		products = dict_detailed_reactions[name_reaction][1]
		products = [x for x in products if x != []] #remove empty lists

		nb_products = len(products)
		cpt = 1
		for p in  products:
			if cpt < nb_products:
				one_reaction += "".join(p) + " + "
			else:
				one_reaction += "".join(p)
			cpt += 1

	return one_reaction

#************************************************

def in_and_out_metabolites(DICT_reactions):

	dict_in_M = {}
	dict_out_M = {}

	for reaction, data in DICT_reactions.items():

		if len(data) == 2: #not reversible

			reactants = data[0]
			reactants = [x for x in reactants if x != []] #remove empty lists

			products = data[1]
			products = [x for x in products if x != []] #remove empty lists

			in_M = re.findall(r'M_.*_e', str(reactants).replace("'","")) #.replace to remove weird apostrophes
			out_M = re.findall(r'M_.*_e', str(products).replace("'",""))

			if len(in_M)!=0:
				dict_in_M[reaction] = in_M
			if len(out_M)!=0:
				dict_out_M[reaction] = out_M

	return dict_in_M, dict_out_M

#************************************************

def save_notReversibleReactions_file(d, dm, dr, drdetailed, output):

	f = open(output,"w")

	for metabolites,reactions in d.items():

		# print("reactions :" ,reactions, "longueur :", len(reactions))

		if len(reactions)==2: #only not reversible reactions

			#be careful with list of metabolites here !!!		
			split_metabolites = re.split(", ",metabolites)

			for m in split_metabolites:

				f.write("\n--------- "+dm[m]+" ("+m+") ----------\n")
				f.write("from :\n")

				for r in reactions[0]:
					f.write("* "+"".join(dr[r])+" ("+r+")\n")
					f.write("\t")

					reaction = pretty_ONE_reaction(dr,drdetailed,r)
					f.write("\n\t"+reaction+"\n")


				f.write("\nto :\n")

				for r in reactions[1]:
					f.write("* "+"".join(dr[r])+" ("+r+")\n")
					f.write("\t")
					reaction = pretty_ONE_reaction(dr,drdetailed,r)
					f.write("\n\t"+reaction+"\n")

	f.close()

#************************************************

def save_commonMetabolites_file(drdetailed1, drdetailed2, output):

	f = open(output,"w")

	list_extr_metabolites_bacteria1 = []
	list_extr_metabolites_bacteria2 = []

	for reaction, extr_metabolites in drdetailed1.items():
		
		if len(extr_metabolites) == 1:
			if extr_metabolites!=[]:
				list_extr_metabolites_bacteria1.extend(extr_metabolites)
		elif  len(extr_metabolites) > 1 :
			if extr_metabolites[0]!=[]:
				list_extr_metabolites_bacteria1.extend(extr_metabolites[0][0])
			if extr_metabolites[1]!=[]:
				list_extr_metabolites_bacteria1.extend(extr_metabolites[1][0])

	for reaction, extr_metabolites in drdetailed2.items():
		if len(extr_metabolites) == 1:
			if extr_metabolites!=[]:
				list_extr_metabolites_bacteria2.extend(extr_metabolites)
		elif  len(extr_metabolites) > 1 :
			if extr_metabolites[0]!=[]:
				list_extr_metabolites_bacteria2.extend(extr_metabolites[0][0])
			if extr_metabolites[1]!=[]:
				list_extr_metabolites_bacteria2.extend(extr_metabolites[1][0])

	common = sorted(set(list_extr_metabolites_bacteria1).intersection(set(list_extr_metabolites_bacteria2)))

	for m in common:
		f.write(m+"\n")

	f.close()


#************************************************
#************** MAIN ****************************
#************************************************

# Load human data
with open('../Database/Recon3D/Homo_Sapiens_reactions_database', 'rb') as f:
	mon_depickler = pickle.Unpickler(f)
	human_reactions = mon_depickler.load()

with open('../Database/Recon3D/Homo_Sapiens_metabolites_database', 'rb') as f:
	mon_depickler = pickle.Unpickler(f)
	human_metabolites = mon_depickler.load()

# dictionnaire : id métabolites -> nom entier métabolites
metabo_bacterias = make_dict_metabolites(sys.argv[1],human_metabolites)
#dictionnaire : id réaction -> nom entier des réactions
reactions_bacterias = make_dict_reactions(sys.argv[1],human_reactions)

#réactions format traditionnel : réactant1 + réactant2 + ... -> produit1 + produit2 + ...
#avec un composé uniquement des id
#l'autre avec les noms entiers

dict_reactions_bacteria1, dict_detailed_reactions_bacteria1 = get_dict_reactions(sys.argv[1],metabo_bacterias)

dict_reactions_bacteria2, dict_detailed_reactions_bacteria2 = get_dict_reactions_homo_sapiens(sys.argv[2],metabo_bacterias,reactions_bacterias)


#merge two dictionaries
dict_detailed_reactions = {}
for key in set().union(dict_detailed_reactions_bacteria1, dict_detailed_reactions_bacteria2):
    if key in dict_detailed_reactions_bacteria1: dict_detailed_reactions.setdefault(key, []).extend(dict_detailed_reactions_bacteria1[key])
    if key in dict_detailed_reactions_bacteria2: dict_detailed_reactions.setdefault(key, []).extend(dict_detailed_reactions_bacteria2[key])

#dictionnaire : id réaction : id métabolites entrants et sortants
dict_in_bacteria1,dict_out_bacteria1= in_and_out_metabolites(dict_reactions_bacteria1)
dict_in_bacteria2,dict_out_bacteria2= in_and_out_metabolites(dict_reactions_bacteria2)

# print(dict_in_bacteria1,dict_in_bacteria2)
# print(dict_out_bacteria1,dict_out_bacteria2)

#listes des métabolites va et vient en
#de bactérie 1 à 2
bacteria1_to_bacteria2 = set(sum(dict_out_bacteria1.values(), [])).intersection(set(sum(dict_in_bacteria2.values(), [])))
#de bactérie 2 à 1
bacteria2_to_bacteria1 = set(sum(dict_out_bacteria2.values(), [])).intersection(set(sum(dict_in_bacteria1.values(), [] )))

#maintenant les associer à la réaction chimique correspondante depuis la bactérie 1 vers la 2 et vice verça
out_R = []
in_R = []
dict_bacteria1_to_bacteria2 = {}
for m in bacteria1_to_bacteria2:
	for u,v in dict_out_bacteria1.items(): 
		if m in v:
			out_R.append("".join(u))
	for u,v in dict_in_bacteria2.items(): 
		if m in v:
			in_R.append("".join(u))
	dict_bacteria1_to_bacteria2[m] = out_R , in_R
	out_R = []
	in_R = []

out_R = []
in_R = []
dict_bacteria2_to_bacteria1 = {}
for m in bacteria2_to_bacteria1:
	for u,v in dict_out_bacteria2.items(): 
		if m in v:
			out_R.append("".join(u))
	for u,v in dict_in_bacteria1.items(): 
		if m in v:
			in_R.append("".join(u))
	dict_bacteria2_to_bacteria1[m] = [out_R, in_R]
	out_R = []
	in_R = []

bacteria1_name = re.split(".xml",re.split("/",sys.argv[1])[-1])[0]
bacteria2_name = re.split(".xml",re.split("/",sys.argv[2])[-1])[0]

if bacteria1_to_bacteria2 != set() and bacteria2_to_bacteria1 != set(): 

	out_R = []
	in_R = []
	dict_bacteria1_to_bacteria2 = {}
	for m in bacteria1_to_bacteria2:
		for u,v in dict_out_bacteria1.items(): 
			if m in v:
				out_R.append("".join(u))
		for u,v in dict_in_bacteria2.items(): 
			if m in v:
				in_R.append("".join(u))
		dict_bacteria1_to_bacteria2[m] = out_R , in_R
		out_R = []
		in_R = []

	out_R = []
	in_R = []
	dict_bacteria2_to_bacteria1 = {}
	for m in bacteria2_to_bacteria1:
		for u,v in dict_out_bacteria2.items(): 
			if m in v:
				out_R.append("".join(u))
		for u,v in dict_in_bacteria1.items(): 
			if m in v:
				in_R.append("".join(u))
		dict_bacteria2_to_bacteria1[m] = [out_R, in_R]
		out_R = []
		in_R = []

	output1 = sys.argv[3]+"/"+bacteria1_name+"_to_"+ bacteria2_name+".txt"
	output2 = sys.argv[3]+"/"+bacteria2_name+"_to_"+ bacteria1_name+".txt"

	# print("ICI",dict_detailed_reactions_bacteria1["Postulated transport reaction"])
	# print("LA",dict_detailed_reactions_bacteria2["Postulated transport reaction"])

	save_notReversibleReactions_file(dict_bacteria1_to_bacteria2, metabo_bacterias, reactions_bacterias,dict_detailed_reactions,output1)

	save_notReversibleReactions_file(dict_bacteria2_to_bacteria1, metabo_bacterias, reactions_bacterias,dict_detailed_reactions,output2)

###Reversible reactions

output = sys.argv[3]+"/"+bacteria1_name+"_"+ bacteria2_name+"_common_extracellular_metabolites.txt"

save_commonMetabolites_file(dict_detailed_reactions_bacteria1,dict_detailed_reactions_bacteria2, output)

"""
	le délire c'est de trouver les sets communs de extracellulaire uniquement en fait. ouaiiis. 
	puis je fais des listes.
	bacteria1_bacteria2_metabolites_extracellular.txt
""" 


