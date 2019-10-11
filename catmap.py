#!/usr/bin/env python

import numpy as np
import pickle
from copy import deepcopy
import os
import sys
from shared_functions import carbon_number, carbon_fragments, inverse_structure

def catmap():
	'writes the reaction in catmap format'
	'load data from already generated pickle files containing reactants, ts, and final products'
	'it also generates a barrier dictionary'
	#barriers = {}
	#barriers['O'], barriers['H'], barriers['C'],barriers['assisted'] = {},{},{},{}
	#rxn_energy = {}
	reactant = sys.argv[1]
	#energies = pickle.load((open(reactant+'-energies.pkl','rb'))) #contains all formation energies
	f = open(reactant+'-catmap.log','w')
	data = pickle.load((open(reactant+'-ts-dict-reduced.pkl','rb')))
	pro_list = data.keys()
	for reactant in data:
		#print reactant
		for fs in data[reactant]:
		      'make sure elements are in products list'
		      if fs[0] not in pro_list and carbon_number(fs[0])[0] != 1:
				print 'element not in product list: ', fs[0]
				continue
		      elif fs[1] not in pro_list and carbon_number(fs[0])[0] != 1 and fs[1] not in ['H','O','OH','O ','O  ','O   ','OH ','OH ','OH   ','OH    ']:
				print '??', fs[1]
                                continue
			
		      else:

			'write each elementary reaction'
			if fs[1] in ['O','O ','O  ','O   ']:
                            if data[reactant][fs].replace('-','') in data:
                                FS = data[reactant][fs].replace('-','')
			        f.write("'"+reactant+'_s + O_s <-> '+data[reactant][fs]+'_s + *_s <-> '+FS+'_s + *_s'+"',"' \n')
				#barriers['O'][data[reactant][fs]] = round(energies[data[reactant][fs]+'_211'] - energies[reactant+'_211'] - energies['O_211'],3)
				#rxn_energy[data[reactant][fs]] = energies[FS+'_211'] - energies['O_211'] - energies[reactant+'_211']
			
				
                            elif inverse_structure(data[reactant][fs].replace('-','')) in data:
                                FS = inverse_structure(data[reactant][fs].replace('-',''))
			        f.write("'"+reactant+'_s + O_s <-> '+data[reactant][fs]+'_s + *_s <-> '+FS+'_s + *_s'+"',"' \n')
				#barriers['O'][data[reactant][fs]] = round(energies[data[reactant][fs]+'_211'] - energies[reactant+'_211'] - energies['O_211'],3)
				#rxn_energy[data[reactant][fs]] = energies[FS+'_211'] - energies['O_211'] - energies[reactant+'_211']
                            else:
                                continue

			elif fs[1] in ['OH','OH ','OH  ','OH   ']:
                            if data[reactant][fs].replace('-','') in data:
                                FS = data[reactant][fs].replace('-','')
			        f.write("'"+reactant+'_s + OH_s <-> '+data[reactant][fs]+'_s + *_s <-> '+FS+'_s + *_s'+"',"' \n')
				#barriers['O'][data[reactant][fs]] = round(energies[data[reactant][fs]+'_211'] - energies[reactant+'_211'] - energies['OH_211'],3)
				#rxn_energy[data[reactant][fs]] = energies[FS+'_211'] - energies['OH_211'] - energies[reactant+'_211']
                            elif inverse_structure(data[reactant][fs].replace('-','')) in data:
                                FS = inverse_structure(data[reactant][fs].replace('-',''))
			        f.write("'"+reactant+'_s + OH_s <-> '+data[reactant][fs]+'_s + *_s <-> '+FS+'_s + *_s'+"',"' \n')
				#barriers['O'][data[reactant][fs]] = round(energies[data[reactant][fs]+'_211'] - energies[reactant+'_211'] - energies['OH_211'],3)
				#rxn_energy[data[reactant][fs]] = energies[FS+'_211'] - energies['OH_211'] - energies[reactant+'_211']
                            else:
                                continue
                        else:
			    'writes elementary reaction for C-H and C-C'
			    f.write("'"+reactant+'_s + *_s <-> '+data[reactant][fs]+'_s + *_s <-> '+fs[0]+'_s + '+fs[1]+'_s'+"',"' \n')
			    if '-C' in data[reactant][fs]:
				    #if data[reactant][fs]+'_211' in energies:
				    	#barriers['C'][data[reactant][fs]] = round(energies[data[reactant][fs]+'_211'] - energies[reactant+'_211'],3)
			    	f.write("'"+reactant+'_s + *_s <-> '+data[reactant][fs]+'_s + *_s <-> '+fs[0]+'_s + '+fs[1]+'_s'+"',"' \n')
				#elif inverse_structure(data[reactant][fs])+'_211' in energies:
				    	#barriers['C'][data[reactant][fs]] = round(energies[inverse_structure(data[reactant][fs])+'_211'] - energies[reactant+'_211'],3)
			    	f.write("'"+reactant+'_s + *_s <-> '+inverse_structure(data[reactant][fs])+'_s + *_s <-> '+fs[0]+'_s + '+fs[1]+'_s'+"',"' \n')
				#else: 
					#print 'ts not in ts dict! ', data[reactant][fs] 
				    	#barriers['C'][data[reactant][fs]] = round(energies[data[reactant][fs]+'_211'] - energies[reactant+'_211'],3)
			    else:
				#if data[reactant][fs]+'_211' in energies:
				    #barriers['H'][data[reactant][fs]] = round(energies[data[reactant][fs]+'_211'] - energies[reactant+'_211'],3)
			    	    f.write("'"+reactant+'_s + *_s <-> '+data[reactant][fs]+'_s + *_s <-> '+fs[0]+'_s + '+fs[1]+'_s'+"',"' \n')
				#elif inverse_structure(data[reactant][fs])+'_211' in energies:
				    #barriers['H'][data[reactant][fs]] = round(energies[inverse_structure(data[reactant][fs])+'_211'] - energies[reactant+'_211'],3)
			    	    f.write("'"+reactant+'_s + *_s <-> '+inverse_structure(data[reactant][fs])+'_s + *_s <-> '+fs[0]+'_s + '+fs[1]+'_s'+"',"' \n')
				#else:
					#print '>>>not in ts dict!', data[reactant][fs]

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
				
				#if TS_O+'_211' and TS_OH+'_211' in energies:
					#continue
					#barriers['assisted'][TS_O] = round(energies[TS_O+'_211'] - energies[reactant+'_211'] - energies['O_211'],3)
					#barriers['assisted'][TS_OH] = round(energies[TS_OH+'_211'] - energies[reactant+'_211'] - energies['OH_211'],3)
				#else:
					#print '>>>> ', TS_O
	#pickle.dump(barriers, open(sys.argv[1]+'-barriers.pkl','wb'),0)
	#pickle.dump(rxn_energy, open(sys.argv[1]+'-O-rxn-energies.pkl','wb'),0)
	#print barriers
	f.close()

catmap()
