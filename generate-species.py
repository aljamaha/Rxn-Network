from copy import deepcopy
import pickle
import sys
from shared_functions import duplicate_products, carbon_number, carbon_bonds

'Objective:'
'generate a dictionary for initial, transition, and final states for each elementary step/ts/fs'

'Output:'
'Returns a dictionary of reactant, transition, and final states. The content is as follow:'
'for C-H/C-C:  [reactant][fs][ts]'
'for C-O/C-OH: [irrelevant][ts][is]'

def OOH_products(specie):
	'function to add O/OH for each carbon that has < 4 bonds'
        'it also creates a dictionary of 3 levels: specie to which attach O/OH >> ts >> is'
        'e.g. ts[CH][CH-O][CH+O]' 
	'this script can be much shorter'
	C_index = carbon_number(specie)[1] #index of each carbon atom
	products = [] #to collect products into a list
	fs_list = []
	ts_dict = {}
	for i,index in enumerate(C_index):
		O,OH = list(specie), list(specie) #split specie name from string to a list so it could easily be adjusted
		if index == max(C_index):
			'modify last carbon'
			location = 'end'
			if carbon_bonds(specie, index, location) == 4:
				continue
			else:
				O.insert(len(specie),'O')
				OH.insert(len(specie),'OH')

				TS = deepcopy(list(O))
				TS.insert(len(specie),'-')
				fs_list.append(["".join(specie),'O'])
				ts_dict["".join(specie),'O'] = "".join(TS)

				TS = deepcopy(list(OH))
				TS.insert(len(specie),'-')
				fs_list.append(["".join(specie),'OH'])
				ts_dict["".join(specie),'OH'] = "".join(TS)

		elif index == min(C_index):
			'modify starting carbon'
			location = 'beginning'
			if carbon_bonds(specie, index, location) == 4:
				continue
			else:
				O.insert(C_index[i+1],'O')
				OH.insert(C_index[i+1],'OH')

				TS = deepcopy(list(O))
				TS.insert(C_index[i+1],'-')
				fs_list.append(["".join(specie),'O'])
				ts_dict["".join(specie),'O '] = "".join(TS)

				TS = deepcopy(list(OH))
				TS.insert(C_index[i+1],'-')
				fs_list.append(["".join(specie),'OH'])
				ts_dict["".join(specie),'OH '] = "".join(TS)
		else:
			'modify middle carbon'
			location = 'middle'
			if carbon_bonds(specie, index, location) == 4:
				continue
			else:
				O.insert(C_index[i+1],'O')
				OH.insert(C_index[i+1],'OH')

				TS = deepcopy(list(O))
				TS.insert(C_index[i+1],'-')
				fs_list.append(["".join(specie),'O'])
				ts_dict["".join(specie),'O   '] = "".join(TS)

				TS = deepcopy(list(OH))
				TS.insert(C_index[i+1],'-')
				fs_list.append(["".join(specie),'OH'])
				ts_dict["".join(specie),'OH   '] = "".join(TS)

		products.append(''.join(O))
		products.append(''.join(OH))

	return products, ts_dict

def CH_products(specie):
        'generates all possible C-C and C-H products from a certain specie'
        specie = list(specie) #convert specie into a list so it I can easily modify the strings
        products_list = []
	fs_list = []
	ts_dict = {}
	for index,elements in enumerate(specie):
                'loop over every element in the reactant speice'
		if elements == 'H':
			'remove a hydrogen atom to create a potential product specie'
			if index < len(specie) - 1:
				#if it is not the last element in the specie
				if specie[int(index)+1] == '3':
					'if it an H3, make it an H2'
					product = deepcopy(specie)
					product[index+1] = '2'
					products_list.append("".join(product))
					fs_list.append(["".join(product),'H'])
					TS = deepcopy(specie)
					TS[index+1] = '2'
					TS.insert(index+2,'-')
					TS.insert(index+3,'H')
					ts_dict["".join(product),'H'] = "".join(TS)
				elif specie[int(index)+1] == '2':
					'if it an H2, make it an H'
					product = deepcopy(specie)
					product[index+1] = ''
					products_list.append("".join(product))
					fs_list.append(["".join(product),'H'])
					TS = deepcopy(specie)
					TS[index+1] = '-'
					TS.insert(index+2,'H')
					ts_dict["".join(product),'H'] = "".join(TS)

                                else:
                                        'remove the H'
                                        product = deepcopy(specie)
                                        product[index] = ''
                                        products_list.append("".join(product))
					fs_list.append(["".join(product),'H'])
					TS = deepcopy(specie)
					TS[index] = '-'
					TS.insert(index+1,'H')
					ts_dict["".join(product),'H'] = "".join(TS)
			else:
				'if it an H, remove it'
				product = deepcopy(specie)
				product[index] = ''
				products_list.append("".join(product))
				fs_list.append(["".join(product),'H'])
				TS = deepcopy(specie)
				TS[index] = '-'
				TS.insert(index+1,'H')
				ts_dict["".join(product),'H'] = "".join(TS)
		elif elements == 'C':
                        'products from breaking the C-C bond'
                        for i in range(index,len(specie)):
                            if specie[i] == 'C':
                                if i == 0:
                                    continue #'if C is first atom, avoids writing the same specie (no C-C cracking)'
                                else:
                                    product = deepcopy(specie)
                                    products_list.append("".join(specie[0:i]))
                                    products_list.append("".join(specie[i:len(specie)]))
				    fs_list.append([("".join(specie[0:i])),("".join(specie[i:len(specie)]))])
				    TS = deepcopy(specie)
				    TS = "".join(specie[0:i])+'-'+"".join(specie[i:len(specie)])
				    ts_dict[("".join(specie[0:i])),("".join(specie[i:len(specie)]))] =  TS
        return products_list, fs_list,ts_dict

fs,ts = {},{}

def all_products(reactant):
	'combines all previous functions to predict all possible products'
	products = [] #storing all the products
	CH =  CH_products(reactant)[0] #generate CH products
	fs[reactant] = CH_products(reactant)[1]
	products = duplicate_products(CH) #remove duplicates from CH products
	products.extend(OOH_products(reactant)[0]) #add O/OH products
	a = CH_products(reactant)[2]
	b = OOH_products(reactant)[1]
	a.update(b) #joining the two dictionaries
	ts[reactant] = a

	return products, fs, ts

def final(reactant):
	'generate a list of all products and a dict of reacants and their products'
	pro_list,pro_dict = [] , {} #compiled list of products and a dict of reactants and products
	pro_list.extend(all_products(reactant)[0]) #creating initial products from there
	pro_dict[reactant] = pro_list #initial creation of the dict of reactants and associated products
	for specie in pro_list:
		if specie in pro_dict:
			continue
		else:
			pro_dict[specie] = all_products(specie)[0]
			pro_list.extend(pro_dict[specie])
	#pickle.dump(pro_dict, open(reactant+"-dictionary.pkl", "wb"), 0) #dump results in a pickle file
	pickle.dump(pro_list, open(reactant+'-products-list.pkl','wb'),0)
	transition_states = all_products(reactant)[2]
	pickle.dump(transition_states,open(reactant+"-ts-dict-full.pkl",'wb'),0)
	return pro_list, pro_dict, transition_states

'lets make this work!'
reactant = sys.argv[1]
p, d, s = final(reactant)
