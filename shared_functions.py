from copy import deepcopy
import pickle
import sys

'these functions are shared by most scripts used in generating the reaction network'

def duplicate_products(products):
    'remove duplicate products, such as double counted CH2'
    all_products, new_list = deepcopy(products),[]
    for product in products:
        if product in new_list:
            continue #'avoid writing the same product twice'
        elif product == '':
            continue #'avoid writing an empty specie'
        else:
            new_list.append(product)
    return new_list

def carbon_number(specie):
    'find the numbe of carbon atoms in the specie and their indicies'
    n , C_index = 0, []
    for index, element in enumerate(specie):
        if element == 'C':
	    C_index.append(index)
            n += 1
    return n, C_index    

def carbon_bonds(specie,carbon_index,location):
	'check the number of bonds C have'
        #if #of bonds = 4, then O/OH should not be added'
	if location == 'beginning':
		bonds = 1 #counts connecting to next carbon
		next_carbon_index = specie[carbon_index+1:len(specie)+1].find('C') + 1 #from the atom next to C to end
	elif location == 'middle':
		bonds = 2 #counts one bond for carbon before and carbon after
		next_carbon_index = specie[carbon_index+1:len(specie)+1].find('C') #from atom next to C to end
		next_carbon_index = next_carbon_index + carbon_index + 1 #bcs search above doesn't start from beg of spec.
	elif location == 'end':
		bonds = 1 #counts connecting to next carbon
		next_carbon_index = len(specie) #when we do a loop, it will stop at the last index

        if next_carbon_index == carbon_index+1: #if it is carbon bonded to another carbon, then we have only one bond
            return bonds
        for i in range(carbon_index,next_carbon_index):
		if i == len(specie) - 1:
			#when it is the last specie
			if specie[i-1:i+1] == 'OH' or specie[i] in ['2','3']:
				# do not double count if it is OH or H2/H3
				continue
			else:
				bonds += 1
		else:
			#if it is not the last specie
			if specie[i:i+2] == 'OH':
        	            bonds +=1
 	                if specie[i] == 'H':
				if specie[i-1] == 'O':
					#do not double count H if we have OH
					continue
		        	elif specie[i+1]=='2':
					bonds +=2
				elif specie[i+1] == '3':
					bonds +=3
				elif specie[i+1] in ['C','O','H']:
					bonds +=1
			elif specie[i] == 'O':
					if specie[i+1] != 'H':
						# it is not OH
						bonds +=1
	return bonds

def carbon_fragments(specie):
   'separate a molecule with more than one carbon atom into individual items (eg. CCHCH2 = C+CH+CH2'
   'output is in the order of the original species'
   n, C_index = carbon_number(specie) #get the number of carbon atoms
   fragment = range(n) #create a list to store fragments
   if carbon_number(specie)[0] == 3:
	for i,index in enumerate(C_index):
		if index == max(C_index):
			'if it is the last carbon'
			fragment[2] = specie[index:len(specie)+1]
		elif index == min(C_index):
			'if it is the first carbon'
			fragment[0] = specie[0:C_index[1]]
		else:
			'if it an intermediete carbon'
			fragment[1] = specie[index:C_index[2]]
   elif carbon_number(specie)[0] == 2:
	for i,index in enumerate(C_index):
		if index == max(C_index):
			fragment[1] = specie[index:len(specie)+1]
		elif index == min(C_index):
			fragment[0] = specie[0:C_index[1]]
   elif carbon_number(specie)[0] == 1:
	return specie
   return fragment

def C_inverse(item):
        'write the inverse structure e.g. CHC becomes CCH'
        if carbon_number(item)[0] == 2:
            new_specie = carbon_fragments(item)[1] + carbon_fragments(item)[0]
        elif carbon_number(item)[0] == 3:
            new_specie = carbon_fragments(item)[2] + carbon_fragments(item)[1] + carbon_fragments(item)[0]
        else:
            new_specie = item
        return new_specie

def convert_structure_to_invserse(item):
	'feed it a carbon specie, and it generates the inverse structure'
	if carbon_number(item)[0] in [0,1]:
		new_specie = item
	elif carbon_number(item)[0] == 2:
		new_specie = carbon_fragments(item)[1] + carbon_fragments(item)[0]
	elif carbon_number(item)[0] == 3:
		new_specie = carbon_fragments(item)[2] + carbon_fragments(item)[1] + carbon_fragments(item)[0]
	return new_specie

def inverse_structure(item):
        'write the inverse structure e.g. CHC becomes CCH'
        if carbon_number(item)[0] == 2:
	  if '-' in item:
		f1 = deepcopy(carbon_fragments(item)[1])
		f2 = deepcopy(carbon_fragments(item)[0])
		f1 = f1.replace('-','')
		f2 = f2.replace('-','')
		new_specie = f1+'-'+f2
		
	  else:
            new_specie = carbon_fragments(item)[1] + carbon_fragments(item)[0]
        elif carbon_number(item)[0] == 3:
            new_specie = carbon_fragments(item)[2] + carbon_fragments(item)[1] + carbon_fragments(item)[0]
        else:
            new_specie = item
        return new_specie
