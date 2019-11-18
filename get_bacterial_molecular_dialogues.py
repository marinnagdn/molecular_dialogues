#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#to get cytosol -> extracellular and extracellular -> cytosol metabolites from a xml metabolism bacteria file

# python3 extract_c_e_metabolites.py /home/marinna/M2/Bioinfo_appliquée_2/embl_gems-master/models/a/acidobacterium/toy_Acidobacterium_ailaaui_PMMR2.xml

import sys
import re

#*************************************

def make_dict_metabolites(xml_file1,xml_file2):

	dict_metabolites = {}

	nameHandle = open(xml_file1,"r")

	for line in nameHandle :

		if re.search("<species.*name=", line):
			split_line = re.split('"',line)
			id = split_line[1]
			name = split_line[3]
			dict_metabolites [id] = name

	nameHandle.close()

	nameHandle = open(xml_file2,"r")

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

def make_dict_reactions(xml_file1, xml_file2):

	dict_reactions = {}

	nameHandle = open(xml_file1,"r")

	for line in nameHandle :

		if re.search("<reaction.*id=", line):
			split_line = re.split('"',line)
			id = split_line[1]
			name = split_line[3]
			dict_reactions[id] = name

	nameHandle.close()

	nameHandle = open(xml_file2,"r")

	for line in nameHandle :

		if re.search("<reaction.*id=", line):
			split_line = re.split('"',line)
			id = split_line[1]
			if id not in dict_reactions.keys():
				name = split_line[3]
				dict_reactions [id] = name

	nameHandle.close()

	return dict_reactions


#*************************************

def get_dict_reactions(xml_file, metabo_bacteria):

	has_found_reaction = False
	has_found_reactant = False
	has_found_product = False
	description_reaction = []
	description_reactants = []
	description_products = []

	dict_reactants_extracellular_M = {}
	dict_reactants_periplasm_M = {}
	dict_reactants_cytosol_M = {}
	dict_products_extracellular_M ={}
	dict_products_periplasm_M = {}
	dict_products_cytosol_M = {}

	dict_reactions = {} #dict with id reactions as keys, and list of lists 1) reactants 2) products as values
	dict_detailed_reactions = {}

	out_M = [] #from cytosol or periplasm to extracellular
	in_M = [] #from extracellular to cytosol or periplasm

	nameHandle = open(xml_file,"r")

	for line in nameHandle :
		if re.search(r"</listOfReactionsQ>",line): #end of reactions
			break
		elif re.search(r"<reaction id=",line): #new reaction
			reaction_id = re.split('"',line)[1]	
			reaction_name = re.split('"',line)[3]	
			# print(reaction_id,reaction_name)
			if not re.search("EX",reaction_id): #delete EX reactions
				has_found_reaction = True #lets get informations of this reaction
		elif has_found_reaction:
			# print (line)
			description_reaction.append(line) #accumule informations of the reactions

		#End of the current reaction : let's analyse it
		if re.search(r"</reaction>",line): 

			for line2 in description_reaction:
				# print(line2)
				#Reactants of the reaction
				if re.search("<listOfReactants>",line2):
					has_found_reactant = True
				elif has_found_reactant:
					description_reactants.append(line2)

				if re.search(r"</listOfReactants>",line2):
					
					for reactant in description_reactants:

						r = re.findall(r"M_.*_p\"",reactant)
						if r != []:
							dict_reactants_periplasm_M["".join(r)[:-1]] = ""


						r = re.findall(r"M_.*_c\"",reactant)
						if r != []:
							dict_reactants_cytosol_M["".join(r)[:-1]] = ""


						r = re.findall(r"M_.*_e\"",reactant)
						if r != []:
							dict_reactants_extracellular_M["".join(r)[:-1]] = ""

					has_found_reactant = False
					description_reactants = []

				#Products of the reaction
				if re.search(r"<listOfProducts>",line2):
					has_found_product = True
				elif has_found_product:
					description_products.append(line2)

				if re.search(r"</listOfProducts>",line2):

					for product in description_products:

						p = re.findall(r"M_.*_p\"",product)
						# p_name = re.split()
						if p != []:
							dict_products_periplasm_M["".join(p)[:-1]] = ""

						p = re.findall(r"M_.*_c\"",product)

						if p != []:
							dict_products_cytosol_M["".join(p)[:-1]] = ""

						p = re.findall(r"M_.*_e\"",product)
						if p != []:
							dict_products_extracellular_M["".join(p)[:-1]] = ""

					has_found_product = False
					description_products = []

			#analysis all metabolite lists.

			#COMING metabolites
			if (len(dict_reactants_extracellular_M.keys()) != 0):
				# in_M.extend(list_reactants_extracellular_M)
				# print(reaction_id)
				dict_reactions[reaction_id] = [list(dict_reactants_extracellular_M.keys()),list(dict_reactants_cytosol_M.keys()), list(dict_reactants_periplasm_M.keys())],[list(dict_products_extracellular_M.keys()),list(dict_products_cytosol_M.keys()),list(dict_products_periplasm_M.keys())]
				dict_detailed_reactions[reaction_name] = [[metabo_bacteria[x] for x in list(dict_reactants_extracellular_M.keys())],  [metabo_bacteria[x] for x in list(dict_reactants_cytosol_M.keys())], [metabo_bacteria[x] for x in list(dict_reactants_periplasm_M.keys())]], [[metabo_bacteria[x] for x in list(dict_products_extracellular_M.keys())],[metabo_bacteria[x] for x in list(dict_products_cytosol_M.keys())],[metabo_bacteria[x] for x in list(dict_products_periplasm_M.keys())]]

			#EXITING metabolites
			if (len(dict_products_extracellular_M.keys()) != 0):
				# out_M.keys().extend(dict_products_extracellular_M.keys())
				dict_reactions[reaction_id] = [list(dict_reactants_extracellular_M.keys()),list(dict_reactants_cytosol_M.keys()), list(dict_reactants_periplasm_M.keys())],[list(dict_products_extracellular_M.keys()),list(dict_products_cytosol_M.keys()),list(dict_products_periplasm_M.keys())]
				dict_detailed_reactions[reaction_name] = [[metabo_bacteria[x] for x in list(dict_reactants_extracellular_M.keys())],  [metabo_bacteria[x] for x in list(dict_reactants_cytosol_M.keys())], [metabo_bacteria[x] for x in list(dict_reactants_periplasm_M.keys())]], [[metabo_bacteria[x] for x in list(dict_products_extracellular_M.keys())],[metabo_bacteria[x] for x in list(dict_products_cytosol_M.keys())],[metabo_bacteria[x] for x in list(dict_products_periplasm_M.keys())]]




			# if len(list_reactants_extracellular_M)!=0 and (len(list_products_cytosol_M)!=0 or len(list_products_periplasm_M)!=0 or len(list_products_extracellular_M)!=0):
				# print("coming")
				# if len(list_products_periplasm_M)!=0:
				# 	in_M.extend(list_products_periplasm_M)
				# if len(list_products_cytosol_M)!=0:
				# 	in_M.extend(list_products_cytosol_M)
				# if len(list_reactants_extracellular_M)!=0:
					# in_M.extend(list_reactants_extracellular_M)

			# elif (len(list_reactants_periplasm_M)!=0 or len(list_reactants_cytosol_M)!=0 or len(list_reactants_extracellular_M)!=0) and len(list_products_extracellular_M)!=0:
				# print("exit")
				# out_M.extend(list_products_extracellular_M)

				# if len(list_reactants_periplasm_M)!=0:
				# if len(list_reactants_cytosol_M)!=0:
				# 	out_M.extend(list_reactants_cytosol_M)

			#empty lists
			dict_reactants_extracellular_M = {}
			dict_reactants_periplasm_M = {}
			dict_reactants_cytosol_M = {}
			dict_products_extracellular_M ={}
			dict_products_periplasm_M = {}
			dict_products_cytosol_M = {}

			description_reaction = [] #erase description
			has_found_reaction = False

	# print("exiting metabolites : ",out_M)
	# print("coming in metabolites : ",in_M)
	# print("list of interest reactions : ", dict_reactions["R_ARGt2r"])

	nameHandle.close()
	# return (out_M, in_M)
	return dict_reactions, dict_detailed_reactions

#*************************************

def pretty_reactions(DICT_reactions):
	for reaction,data in DICT_reactions.items():

		print (reaction,":", end =" ")

		reactants = data[0]
		# print(reactants)
		reactants = [x for x in reactants if x != []] #remove empty lists

		nb_reactants = len(reactants)
		cpt = 1
		for r in  reactants:
			if cpt < nb_reactants:
				print(''.join(r),"+",end =" ")
			else:
				print(''.join(r),end =" ")
			cpt += 1

		print ("->", end =" ")

		products = data[1]
		products = [x for x in products if x != []] #remove empty lists

		nb_products = len(products)
		cpt = 1
		for p in  products:
			if cpt < nb_products:
				print(''.join(p),"+",end =" ")
			else:
				print(''.join(p))
			cpt += 1

#***************************************

def pretty_ONE_reaction(dr, DICT_reactions, id_reaction):

	one_reaction = ""

	name_reaction = dr[id_reaction]

	# print (name_reaction,":", end =" ")

	reactants = DICT_reactions[name_reaction][0]
	# print(reactants)
	reactants = [x for x in reactants if x != []] #remove empty lists

	nb_reactants = len(reactants)
	cpt = 1
	for r in  reactants:
		if cpt < nb_reactants:
			# print(''.join(r),"+",end =" ")
			one_reaction += "".join(r) + " + "
		else:
			# print(''.join(r),end =" ")
			one_reaction += "".join(r)
		cpt += 1

	# print ("->", end =" ")
	one_reaction += " -> "

	products = DICT_reactions[name_reaction][1]
	products = [x for x in products if x != []] #remove empty lists

	nb_products = len(products)
	cpt = 1
	for p in  products:
		if cpt < nb_products:
			# print(''.join(p),"+",end =" ")
			one_reaction += "".join(p) + " + "
		else:
			# print(''.join(p))
			one_reaction += "".join(p)
		cpt += 1

	return one_reaction

#************************************************

def in_and_out_metabolites(DICT_reactions):

	dict_in_M = {}
	dict_out_M = {}


	for reaction, data in DICT_reactions.items():

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

# def pretty_results(d, dm, dr1, dr2, drdetailed1, drdetailed2):

# 	for m,reactions in d.items():

# 		print("\n--------- ",dm[m]," (",m,") ----------\n",sep="")
# 		print("from :")

# 		for r in reactions[0]:
# 			print("* ",dr1[r]," (",r,")\n",sep="")
# 			print("\t",end="")
# 			pretty_ONE_reaction(dr1,drdetailed1,r)


# 		print("\nto :")

# 		for r in reactions[1]:
# 			print("* ",dr2[r]," (",r,")\n",sep="")
# 			print("\t",end="")
# 			pretty_ONE_reaction(dr2,drdetailed2,r)

#************************************************

def save_results_file(d, dm, dr, drdetailed1, drdetailed2, output):

	f = open(output,"w")

	for metabolites,reactions in d.items():

		#be careful with list of metabolites here !!!		
		split_metabolites = re.split(", ",metabolites)

		for m in split_metabolites:

			f.write("\n--------- "+dm[m]+" ("+m+") ----------\n")
			f.write("from :\n")

			for r in reactions[0]:
				f.write("* "+dr[r]+" ("+r+")\n")
				f.write("\t")
				reaction = pretty_ONE_reaction(dr,drdetailed1,r)
				f.write("\n\t"+reaction+"\n")


			f.write("\nto :\n")

			for r in reactions[1]:
				f.write("* "+dr[r]+" ("+r+")\n")
				f.write("\t")
				reaction = pretty_ONE_reaction(dr,drdetailed2,r)
				f.write("\n\t"+reaction+"\n")

	f.close()


#************************************************
#************** MAIN ****************************
#************************************************

#dictionnaire : id métabolites -> nom entier métabolites
metabo_bacterias = make_dict_metabolites(sys.argv[1], sys.argv[2])
#dictionnaire : id réaction -> nom entier des réactions
reactions_bacterias = make_dict_reactions(sys.argv[1],sys.argv[2])

#réactions format traditionnel : réactant1 + réactant2 + ... -> produit1 + produit2 + ...
#avec un composé uniquement des id
#l'autre avec les noms entiers
dict_reactions_bacteria1, dict_detailed_reactions_bacteria1 = get_dict_reactions(sys.argv[1],metabo_bacterias)
dict_reactions_bacteria2, dict_detailed_reactions_bacteria2 = get_dict_reactions(sys.argv[2],metabo_bacterias)

#dictionnaire : id réaction : id métabolites entrants et sortants
dict_in_bacteria1,dict_out_bacteria1= in_and_out_metabolites(dict_reactions_bacteria1)
dict_in_bacteria2,dict_out_bacteria2= in_and_out_metabolites(dict_reactions_bacteria2)

# print(dict_in_bacteria1)
# print(dict_out_bacteria1)
# print(dict_in_bacteria2)
# print(dict_out_bacteria2)

# pretty_reactions(dict_detailed_reactions_bacteria1)

#listes des métabolites va et vient en
#de bactérie 1 à 2
bacteria1_to_bacteria2 = set(sum(dict_out_bacteria1.values(), [])).intersection(set(sum(dict_in_bacteria2.values(), [])))
#de bactérie 2 à 1
bacteria2_to_bacteria1 = set(sum(dict_out_bacteria2.values(), [])).intersection(set(sum(dict_in_bacteria1.values(), [] )))

# print(bacteria1_to_bacteria2)
# print(bacteria2_to_bacteria1)

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

# pretty_results(dict_bacteria1_to_bacteria2, metabo_bacteria1, reactions_bacteria1, reactions_bacteria2,dict_detailed_reactions_bacteria1,dict_detailed_reactions_bacteria2) #l'un ou l'autre metabo_bacteria peu importe comme métabo en commun

# pretty_results(dict_bacteria2_to_bacteria1, metabo_bacteria2, reactions_bacteria2, reactions_bacteria1,dict_detailed_reactions_bacteria2,dict_detailed_reactions_bacteria1) 

#enregistrement des résultats
bacteria1_name = re.split(".xml",re.split("/",sys.argv[1])[-1])[0]
bacteria2_name = re.split(".xml",re.split("/",sys.argv[2])[-1])[0]
# bacteria1_path = re.split(".xml",sys.argv[1])[0]
# bacteria2_path = re.split(".xml",sys.argv[2])[0]

output1 = sys.argv[3]+"/"+bacteria1_name+"_to_"+ bacteria2_name+".txt"
output2 = sys.argv[3]+"/"+bacteria2_name+"_to_"+ bacteria1_name+".txt"

save_results_file(dict_bacteria1_to_bacteria2, metabo_bacterias, reactions_bacterias,dict_detailed_reactions_bacteria1,dict_detailed_reactions_bacteria2, output1)

save_results_file(dict_bacteria2_to_bacteria1, metabo_bacterias, reactions_bacterias,dict_detailed_reactions_bacteria2,dict_detailed_reactions_bacteria1,output2)


