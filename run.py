import os, sys, pickle
from copy import deepcopy
from generate_species import *
from reduce_species import *
from catmap import *
from post_process import *

'Input:'
reactant = sys.argv[1]	#name of the hydrocarbon reactant

'''
Step 1: Generate a dictionary for initial, transition, and final states for each elementary step/ts/fs

Output:
Returns a dictionary of reactant, transition, and final states. The content is as follow:
	for C-H/C-C:  [reactant][fs][ts]
	for C-O/C-OH: [irrelevant][ts][is]
List of all species in the reaction network
'''

print 'generating reactions list ...'
final(reactant)

'''
Step 2: Eliminates repeated products and products with the same strucutre but inversed (CH3CH/CHCH3)

Output:
Returns a dictionary where repeated products are eliminated
'''

print 'eliminating repeated reactions ...'

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
Step 3: wrtie reactions list in CatMap format

Output:
reaction list printed in CatMap input format
'''

print 'writing reactions in catmap format ...'
catmap()

'''
Step 4: remove repeated rxns in CatMap file

Output: 
reaction list printed in CatMap input format without repeated reactions
'''

print 'removing repeated catmap rxns ...'
post_process(reactant)

print('*** Done ***')
