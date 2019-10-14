#!/usr/bin/env python

from ase.atoms import string2symbols
import io, sys, os, pickle, sys
from ase import io
from ase import Atom, Atoms
from ase.atoms import string2symbols
import numpy as np
from copy import deepcopy
from shared_functions import carbon_number, carbon_fragments, inverse_structure

'generates a dictionary of all adsorbates I have and their energies based on existing calculations'

'general items'
wd = '/nfs/slac/g/suncatfs/aljamaha/scaling/C3H6-scaling/Pd/surface-211'
run_dir = os.getcwd()
os.chdir(wd)
cwd = os.getcwd()

'reference energies'
ref_dict = {}
ref_dict['H'] = 0.5*-32.942
ref_dict['O'] = -496.287- 2*ref_dict['H']
ref_dict['C'] = -231.57 - 4*ref_dict['H']
ref_dict['211'] = io.read(wd+'/clean/qn.traj').get_potential_energy()

print ref_dict

def raw_energy_dict(folders):
    'returns a dict with the raw energy of each each species [species in my main folder]'
    energy_dict = {}
    for folder in folders:
	os.chdir(wd)
	if '.' in folder:
		continue
        elif '-' in folder:
            continue
	else:
	   if os.path.exists(folder):
		os.chdir(folder)
		if os.path.exists('qn.traj'):
			atoms = io.read('qn.traj').get_potential_energy()
                        if 'clean' in folder: #for clean surfaces
    			    energy_dict['slab_211'] = atoms
                        elif '_g' in folder: #for gas phase species
    			    energy_dict[folder] = atoms
			else:
    			    energy_dict[folder+'_211'] = atoms
    return energy_dict

'transition state'
def ts_dict(reactant,folders,run_dir):
        'predicts TS energy based on BEP scaling [gives raw energy of TS]'
	'update BEP line slope and intercept'
        ts_raw_energy = {} #dict to store energies in
        os.chdir(run_dir) #dir where we have pickle files
	rxn_dict = pickle.load((open(reactant+'-ts-dict-reduced.pkl','rb'))) #load dict for each rxn
        energies = raw_energy_dict(folders) #loads available raw energies of adsorbates
        rxn_dict.keys() #species available in the ts dictionary [anything else is not part of my rxn network]
        fs_not_in_ts_dict, barriers = [],{}

	for specie in rxn_dict:
		for fs in rxn_dict[specie].keys():
                  'loop over each elementary step in my rxn'
                  if fs[0] not in rxn_dict.keys():
                        'some ts for some reason yield TS not in the dict'
                        'In this elementary step, rxn yields a fs not in my ts dict'
                        if fs[0] not in fs_not_in_ts_dict:
                            fs_not_in_ts_dict.append(fs[0])
                  else:
                    FS1 = fs[0]
                    FS2 = fs[1]
                    FS_string = ''
                    if ' ' in FS2:
                        FS2 = FS2.replace(' ','') #some fs are written as 'O  ', so this renames them to 'O' [removing spaces]
                    'calculate TS based on FS, and then tag them to raw energy of IS'

		    #'C-H bonds'
                    if '-H' in rxn_dict[specie][fs]:
			FS = energies[FS1+'_211'] + energies[FS2+'_211'] 
			composition = string2symbols(specie.replace('-',''))
			for atom in composition:
				FS -= ref_dict[atom]
			FS -= 2*energies['slab_211']
			TS = 0.85*FS + 0.83 
                    	ts_raw_energy[rxn_dict[specie][fs]+'_211'] = round(TS,3)

			'calculate O- and OH-assisted paths'
			FS_O  = energies[FS1+'_211'] + energies['OH_211']
			FS_OH = energies[FS1+'_211'] + energies['H2O_211']
			for atom in composition:
				FS_O  -= ref_dict[atom]
				FS_OH -= ref_dict[atom]
			FS_O  -= ref_dict['O'] + 2*energies['slab_211']
			FS_OH -= ref_dict['O'] + ref_dict['H'] + 2*energies['slab_211']
			TS_O  = 0.56*FS_O + 2.85
			TS_OH = 0.857*FS_OH + 0.87
			index = rxn_dict[specie][fs].find('-')
			TS_O_name = list(rxn_dict[specie][fs])
			TS_O_name.insert(index+2,'--O')
			TS_O_name = "".join(TS_O_name)
			TS_OH_name = list(rxn_dict[specie][fs])
			TS_OH_name.insert(index+2,'--OH')
			TS_OH_name = "".join(TS_OH_name)
                    	ts_raw_energy[TS_O_name+'_211'] = round(TS_O,3)
                    	ts_raw_energy[TS_OH_name+'_211'] = round(TS_OH,3)

		    #'C-C bonds'
		    elif '-C' in rxn_dict[specie][fs]:
			FS = energies[FS1+'_211'] + energies[FS2+'_211']
			composition = string2symbols(specie.replace('-',''))
			for atom in composition:
				FS -= ref_dict[atom]
			FS -= 2*energies['slab_211']
			TS = 0.939*FS + 1.329
                    	ts_raw_energy[rxn_dict[specie][fs]+'_211'] = round(TS,3)
                    else:
                        'calculations for -OH and -O TS'
                        if  rxn_dict[specie][fs].replace('-','') in rxn_dict.keys(): 
                            'IS/FS naming is different in -O/OH. Here I renamed TS by removing - to check if it is in ts dict'
                            FS_string = rxn_dict[specie][fs].replace('-','') #if it is in rxn_dict, I get a new name for FS
                        elif inverse_structure(rxn_dict[specie][fs].replace('-','')) in rxn_dict.keys():
                            'if it is not in ts dict, I check if the inverse is'
                            FS_string = inverse_structure(rxn_dict[specie][fs].replace('-',''))
                        else:
                            if rxn_dict[specie][fs].replace('-','') not in fs_not_in_ts_dict:
                                fs_not_in_ts_dict.append(rxn_dict[specie][fs].replace('-',''))
                        if FS_string == '':
                            'the fact that I do not have a name to FS string means that my FS is not in ts dict'
                            continue

		 	#'C-O bonds'
                        elif '-O' in rxn_dict[specie][fs]:
                            'structure of FS in O/OH is different from C-H/C-C'
                            'two species are combining to create one specie'
                            FS = energies[FS1+'_211'] + energies[FS2+'_211']
			    composition = string2symbols(FS_string.replace('-',''))
			    for atom in composition:
				FS -= ref_dict[atom]
			    FS -= 2*energies['slab_211']
			    TS = 0.969*FS + 1.146
                            ts_raw_energy[rxn_dict[specie][fs]+'_211'] = round(TS,3)
                            'I subtracted the energy of the slab here to make it in a similar nomencluature to C-H/C-C'
        print 'FS is not in my ts dict: ', fs_not_in_ts_dict
        return ts_raw_energy

'convert raw energies to reference energies'
def get_formation_energies(energy_dict,ref_dict):
    formation_energies = {}
    for key in energy_dict.keys(): #iterate through keys
        E0 = energy_dict[key] #raw energy
        name,site = key.split('_') #split key into name/site
        if 'slab' not in name: #do not include empty site energy (0)
	    if 'N' in name:
		'avoiding N2O_g etc.'
		continue
            elif site == '211':
                    E0 -= ref_dict[site] #subtract slab energy if adsorbed
            #remove - from transition-states
            formula = name.replace('-','')
            #get the composition as a list of atomic species
            composition = string2symbols(formula)
            #for each atomic species, subtract off the reference energy
            for atom in composition:
                E0 -= ref_dict[atom]
            #round to 3 decimals since this is the accuracy of DFT
            E0 = round(E0,3)
            formation_energies[key] = E0
            #print key, formation_energies[key]
    return formation_energies

'generate energies.txt'
def make_input_file(file_name,energy_dict,frequency_dict,run_dir):
    os.chdir(run_dir)
    header = '\t'.join(['surface_name','site_name',
                        'species_name','formation_energy',
                        'frequencies','reference'])
    lines = [] #list of lines in the output
    for key in energy_dict.keys(): #iterate through keys
        E = energy_dict[key] #raw energy
        name,site = key.split('_') #split key into name/site
        if 'slab' not in name: #do not include empty site energy (0)
            #frequency = frequency_dict[key]
            frequency = []
            if site == 'gas':
                surface = None
            else:
                surface = 'Pd'
            outline = [surface,site,name,E,frequency,'Input File Tutorial.']
            line = '\t'.join([str(w) for w in outline])
            lines.append(line)

    lines.sort() #The file is easier to read if sorted (optional)
    lines = [header] + lines #add header to top
    input_file = '\n'.join(lines) #Join the lines with a line break

    input = open(file_name,'w') #open the file name in write mode
    input.write(input_file) #write the text
    input.close() #close the file

    print 'Successfully created input file'

'generate a list of energies I have'
os.system("ls > tmp")
energies = [line.rstrip('\n') for line in open('tmp')] #all the species in the main folder
os.system('rm tmp')
os.chdir(cwd)

'lets go!'
reactant = sys.argv[1]
raw_energies   = raw_energy_dict(energies) #generate a dict of adsorption energies
#raw_energies.update(ts_raw_energy) #combine dict of ads and ts energies
energy_dict = get_formation_energies(raw_energies,ref_dict) 
energy_dict.update(ts_dict(reactant,energies,run_dir)) #generate a dict on ts energies
make_input_file(reactant+'-energies.txt',energy_dict,[],run_dir)
pickle.dump(energy_dict, open(reactant+'-energies.pkl','wb'),0)

'large unit cells'
'I think I resolved the issue of how some species were on a large unit cell'
'This, however, indicates if I have somethig wrong'
large = []
for i in energy_dict:
    if np.abs(energy_dict[i])>10:
            name,site = i.split('_')
            large.append(name)
print 'formation energy is wrong (too high) ', large
