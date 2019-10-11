#!/usr/bin/env python

import numpy as np
import pickle
from copy import deepcopy
import os
import sys
from shared_functions import carbon_number, carbon_fragments, inverse_structure

'''
writes the reaction lists in catmap format

Input:
pkl file containing reactants, ts, and final products

Output:
List of elementary steps in CatMap format
'''


def catmap():
	reactant = sys.argv[1]
	f = open(reactant+'-catmap.log','w')
	data = pickle.load((open(reactant+'-ts-dict-reduced.pkl','rb')))
	pro_list = data.keys()
	for reactant in data:
		for fs in data[reactant]:
		      'make sure elements are in products list'
		      if fs[0] not in pro_list and carbon_number(fs[0])[0] != 1:
				#print 'element not in product list: ', fs[0]
				continue
		      elif fs[1] not in pro_list and carbon_number(fs[0])[0] != 1 and fs[1] not in ['H','O','OH','O ','O  ','O   ','OH ','OH ','OH   ','OH    ']:
				#print '??', fs[1]
                                continue
			
		      else:
			'write each elementary reaction'
			if fs[1] in ['O','O ','O  ','O   ']:
                            if data[reactant][fs].replace('-','') in data:
                                FS = data[reactant][fs].replace('-','')
			        f.write("'"+reactant+'_s + O_s <-> '+data[reactant][fs]+'_s + *_s <-> '+FS+'_s + *_s'+"',"' \n')
				
                            elif inverse_structure(data[reactant][fs].replace('-','')) in data:
                                FS = inverse_structure(data[reactant][fs].replace('-',''))
			        f.write("'"+reactant+'_s + O_s <-> '+data[reactant][fs]+'_s + *_s <-> '+FS+'_s + *_s'+"',"' \n')

			elif fs[1] in ['OH','OH ','OH  ','OH   ']:
                            if data[reactant][fs].replace('-','') in data:
                                FS = data[reactant][fs].replace('-','')
			        f.write("'"+reactant+'_s + OH_s <-> '+data[reactant][fs]+'_s + *_s <-> '+FS+'_s + *_s'+"',"' \n')
                            elif inverse_structure(data[reactant][fs].replace('-','')) in data:
                                FS = inverse_structure(data[reactant][fs].replace('-',''))
			        f.write("'"+reactant+'_s + OH_s <-> '+data[reactant][fs]+'_s + *_s <-> '+FS+'_s + *_s'+"',"' \n')
                        else:
			    'writes elementary reaction for C-H and C-C'
			    f.write("'"+reactant+'_s + *_s <-> '+data[reactant][fs]+'_s + *_s <-> '+fs[0]+'_s + '+fs[1]+'_s'+"',"' \n')
			    if '-C' in data[reactant][fs]:
			    		f.write("'"+reactant+'_s + *_s <-> '+data[reactant][fs]+'_s + *_s <-> '+fs[0]+'_s + '+fs[1]+'_s'+"',"' \n')
			    else:
			   	'Dehydrogenation steps'
			    	f.write("'"+reactant+'_s + *_s <-> '+data[reactant][fs]+'_s + *_s <-> '+fs[0]+'_s + '+fs[1]+'_s'+"',"' \n')

			    'writes elementary rxns for C-H [O/OH-assisted]'
			    if fs[1] == 'H':
				index = data[reactant][fs].find('-')
				TS_O =  list(data[reactant][fs])
				TS_O.insert(index+2,'--O')
				TS_O = "".join(TS_O)

				TS_OH =  list(data[reactant][fs])
				TS_OH.insert(index+2,'--OH')
				TS_OH = "".join(TS_OH)

			    	f.write("'"+reactant+'_s + O_s <-> '+TS_O+'_s + *_s <-> '+fs[0]+'_s + OH_s'+"',"' \n')
			    	f.write("'"+reactant+'_s + OH_s <-> '+TS_OH+'_s + *_s <-> '+fs[0]+'_s + H2O_g + *_s'+"',"' \n')
				
	f.close()

catmap()
