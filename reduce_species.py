#!/usr/bin/env python

import pickle
from copy import deepcopy
import os
import sys 
from shared_functions import *


#'Objective:'
#'eliminates repeated products and products with the same strucutre but inversed (CH3CH/CHCH3)'

#reactant = sys.argv[1]

def load_product_list(reactant):
	'load data from already generated pickle files for products list'
	pro_list = pickle.load((open(reactant+'-ts-dict-full.pkl','rb')))
	#print 'original length = ', len(pro_list)
	return pro_list

def inverse_structure(pro_list):
	'eliminate products with same strucutre [e.g. CH2CH and CHCH2]'
	original = deepcopy(pro_list)
	for item in pro_list: #I don't think I need to make it in reversed
		if carbon_number(item)[0] == 2:
			new_specie = carbon_fragments(item)[1] + carbon_fragments(item)[0]
			if new_specie in original:
				if carbon_fragments(item)[1] == carbon_fragments(item)[0]: #e.g. CH2CH2
					continue
				else:
					del original[item]
		elif carbon_number(item)[0] == 3:
			new_specie = carbon_fragments(item)[2] + carbon_fragments(item)[1] + carbon_fragments(item)[0]
			if new_specie in original:
				if carbon_fragments(item)[2] == carbon_fragments(item)[0]:
					continue
				else:
					del original[item]
	return original

def unrealistic_3O_products(pro_list):
	'eliminate product wwhere carbon is bonded to 3 oxygen '
	unrealistic = deepcopy(pro_list)
	for item in pro_list:
		for specie in carbon_fragments(item):
			i = 0
			for element in list(specie):
				if element == 'O':
					i += 1
			if i > 2:
				del unrealistic[item]
				
				break
	#print 'after removing unrealistic products with 3O = ', len(unrealistic)
	return unrealistic

def removing_3O_from_C1(pro_list):
	'removing specie with 3O from a single carbon molecule'
	unrealistic = deepcopy(pro_list)
	for item in pro_list:
		if carbon_number(item)[0] == 1:
			i = 0
			for element in list(item):
				if element == 'O':
					 i += 1
			if i > 2:
				del unrealistic[item]
	#print 'after removing unrealistic 3O products from C1 = ', len(unrealistic)
	return unrealistic

def unlikely_products(pro_list):
	'remove products that are unlikely to form'
	unrealistic = deepcopy(pro_list)
	for item in pro_list:
	   if carbon_number(item)[0] > 1:
		for specie in carbon_fragments(item):
                        if 'COOHH' in specie:
				del unrealistic[item]
				break
			if 'OO' in specie:
				if 'OOH' in specie:
					continue
				else:
					del unrealistic[item]
					break
			if 'CH2O' in specie:
				del unrealistic[item]
				break
			if 'OHOH' in specie:
				del unrealistic[item]
				break
			if 'HOHO' in specie:
				del unrealistic[item]
				break		

	#print 'After removing unlikely items = ', len(unrealistic)
	return unrealistic

def remove_HOOH(pro_list):
	'remove C2/C3 specie with HOOH'
	remove_HOOH = deepcopy(pro_list)
	for item in pro_list:
		if carbon_number(item)[0] > 1:
			for specie in carbon_fragments(item):
				if 'HOOH' in specie:
					del remove_HOOH[item]
					break

	#print 'after removing HOOH = ', len(remove_HOOH)
	return remove_HOOH

def middle_carbon_2O(pro_list):
	'remove middle carbon attached to more than one oxygen'
	cut = deepcopy(pro_list)
	for item in pro_list:
		if carbon_number(item)[0] == 3:
			middle_carbon = item[carbon_number(item)[1][1]:carbon_number(item)[1][2]]
			i = 0
			for element in middle_carbon:
				if element == 'O':
					i +=1
				if i > 1:
					del cut[item]				
					break
	#print 'after removing middle carbon with 2O = ', len(cut)
	return cut

def rename_fs(pro_list):
	'renames fs species so we have no repeated species'
	renamed_list = deepcopy(pro_list)
	for item in pro_list:
		for i, fs in enumerate(pro_list[item].keys()):
			'if a specie in the reactants dict, do nothing'
			if fs[0] in pro_list:
				continue
			else:
				'if a specie not in the reactants dict, rename it'
				new_specie = convert_structure_to_invserse(fs[0])
				renamed_list[item][(new_specie,fs[1])] = renamed_list[item].pop((fs[0],fs[1])) #rename the key in the dict of fs to the one in reactants dict

	'I probably dont need to separate this into two fragments'
	renamed_list_2 = deepcopy(renamed_list)
	for item in renamed_list:
		'do the same for the second fragment in the final state'
		for i, fs in enumerate(pro_list[item].keys()):
			if fs[1] in pro_list:
				continue
			elif fs[1] == 'H':
				continue
			else:
                                #print fs[0],fs[1]
                                if fs[0] or fs[1] in ['COHOHOH','COHOH','COHOO','COHOOH']: #for some reason it is not in products list
                                    continue
                                else:
				    new_specie = convert_structure_to_invserse(fs[1])
				    renamed_list_2[item][(fs[0], new_specie)] = renamed_list_2[item].pop((fs[0],fs[1]))
	return renamed_list_2

def repeated_fs(pro_list):
	'eliminiate two reactions with the same final state'
	new_list = deepcopy(pro_list)
	for item in pro_list:
		for fs in pro_list[item].keys():
			'if two final states are the same, but with positions exchanged'
			new_fs = (fs[1],fs[0])
			if new_fs in pro_list[item].keys():
				if fs[1] == fs[0]:
					continue
				else:
					del pro_list[item][fs]

	'calculate the number of elementary reactions'
	return pro_list

def rename_COHO_to_COOH(pro_list):
	'renames species/ts/fs species from COHO to COOH'

        'renaming fs'
	renamed_list = deepcopy(pro_list)
	for item in pro_list:
	    for i, fs in enumerate(pro_list[item].keys()):
                if carbon_number(fs[0])[0] > 1:
			if 'COHO' in fs[0]:
                            new_specie = fs[0].replace('COHO','COOH')
			    renamed_list[item][(new_specie,fs[1])] = renamed_list[item].pop((fs[0],fs[1])) #rename the key in the dict of fs to the one in reactants dict
	renamed_list_2 = deepcopy(renamed_list)
	for item in renamed_list:
	    'do the same for the second fragment in the final state'
	    for i, fs in enumerate(pro_list[item].keys()):
                if carbon_number(fs[1])[0] > 1:
			if 'COHO' in fs[1]:
                            new_specie = fs[1].replace('COHO','COOH')
			    renamed_list_2[item][(fs[0], new_specie)] = renamed_list_2[item].pop((fs[0],fs[1]))

        'renaming the species (keys in dictionary)'
        renamed_list_3 = deepcopy(renamed_list_2)
        for item in renamed_list_2:
            if 'COHO' in item:
                new_key_name = item.replace('COHO','COOH')
                renamed_list_3[new_key_name] = renamed_list_3.pop(item)

        'renaming TS'
        renamed_list_4 = deepcopy(renamed_list_3)
        for item in renamed_list_3:
            for key in renamed_list_3[item]:
                if 'COHO' in renamed_list_3[item][key]:
                    new = renamed_list_3[item][key].replace('COHO','COOH')
                    renamed_list_4[item][key] = new
        
	#print 'after renaming COHO to COOH ', len(renamed_list_4)

	return renamed_list_4


def CHO_middle_carbon(pro_list):
	'finds and deletes ads with a middle carbon (CHO) where the total # of bonds exceeds the amount the specie can have'
	'I modified it to incluce all species with CHO'
	new_list = deepcopy(pro_list)
	d = 0
	for item in new_list:
		if carbon_number(item)[0] == 3:
			if carbon_fragments(item)[1] == 'CHO':
						del pro_list[item]
	#print 'length after CHO middle carbon =', len(pro_list)
	return pro_list

def CHOH_middle_carbon(pro_list):
	'finds and deletes ads with a middle carbon (CHOH) attached to a carbon bonded to anyhting other than C/CH/COH '
	'this is probably not a great assumption. I need to revisit it later'
	'I made it now so that any CHOH is not going to form, unless it has a terminal C bound to a hollow site'
	'this is due to H being shared to a terminal carbon (if it can accommodate an additional bond. The other occation where structure is fav, then OH pointing out, which I find not very likely'
	new_list = deepcopy(pro_list)
	for item in new_list:
		if carbon_number(item)[0] == 3:
			if carbon_fragments(item)[1] == 'CHOH':
				if carbon_fragments(item)[0] == 'C':
					continue
				elif carbon_fragments(item)[2] == 'C':
					continue
				else:
					del pro_list[item]
	#print 'length after CHOH middle carbon =', len(pro_list)
	return pro_list

def CO_middle_carbon(pro_list):
	'finds and deletes ads with a middle carbon (CO) attached to a carbon bonded to anyhting other than C/CH/COH '
	'this is probably not a great assumption. I need to revisit it later'
	new_list = deepcopy(pro_list)
	for item in new_list:
		if carbon_number(item)[0] == 3:
			if carbon_fragments(item)[1] == 'CO':
				if carbon_fragments(item)[0] == 'C' or carbon_fragments(item)[0] == 'CH' or carbon_fragments(item)[0] == 'COH':
					continue
				elif carbon_fragments(item)[2] == 'C' or carbon_fragments(item)[2] == 'CH' or carbon_fragments(item)[2] == 'COH':
					continue
				else:
					del pro_list[item]
	#print 'length after CO middle carbon =', len(pro_list)
	return pro_list

def CHO_tail_carbon(pro_list):
	'finds CHO that is at the tail connected to a carbon that is attached to anything'
	new_list = deepcopy(pro_list)
	for item in new_list:
		if carbon_number(item)[0] == 3:
			if carbon_fragments(item)[0] == 'CHO' or  carbon_fragments(item)[2] == 'CHO':
				if carbon_fragments(item) == 'C':
					continue
				else:
					del pro_list[item]
	#print 'length after CHO at tail connected to a C-X =', len(pro_list)
	return pro_list

def check_energy(pro_list):
	'check if I have the energy of all species'
	os.chdir('/nfs/slac/g/suncatfs/aljamaha/scaling/C3H6-scaling/Pd/surface-211')
	os.system("ls > tmp")
	folders = [line.rstrip('\n') for line in open('tmp')] #adsorbates I have energy calculations
	List_of_ads = []
	total = 0
	for ads in pro_list:
		'how many energies do I have?'
		if ads in folders:
			total += 1
                elif C_inverse(ads) in folders:
                    total += 1
                #else:
		     #print ads
	#print 'total adsorbates I have energy calculations for = ', total
        #print 'missing calc = ', len(pro_list) - total

def adjacent_C_with_O(pro_list):
	new_list = deepcopy(pro_list)
	for item in new_list:
		if carbon_number(item)[0] == 2:
			if 'O' in carbon_fragments(item)[0] and 'O' in carbon_fragments(item)[1]:
				del pro_list[item]
		elif carbon_number(item)[0] == 3:			
			if 'O' in carbon_fragments(item)[0] and 'O' in carbon_fragments(item)[1]:
				del pro_list[item]
			elif 'O' in carbon_fragments(item)[1] and 'O' in carbon_fragments(item)[2]:
				del pro_list[item]

	#print 'length after adjacent oxygens are removed', len(pro_list)
	return pro_list

'''
loaded_product_list = load_product_list(reactant)
a = inverse_structure(loaded_product_list)
a = rename_COHO_to_COOH(a)
a = rename_fs(a)
a = repeated_fs(a)
a = removing_3O_from_C1(a)
a = unrealistic_3O_products(a)
a = remove_HOOH(a)
a = middle_carbon_2O(a)
a = unlikely_products(a)
a = CHO_middle_carbon(a)
a = CHOH_middle_carbon(a)
a = CO_middle_carbon(a) 
a = CHO_tail_carbon(a)
a = adjacent_C_with_O(a)

pickle.dump(a,open(reactant+"-ts-dict-reduced.pkl",'wb'),0)
'''
