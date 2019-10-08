import os
import sys

reactant = sys.argv[1]
print '='*20, 'generate species'
os.system('python generate-species.py '+reactant) #dict of is/ts/fs
print '='*20, 'reduce species'
os.system('python reduce-species.py '+reactant) #eliminate redundant speciesa
print 'does this also include removing O in adjacent carbons?'
print '='*20, 'catmap'
os.system('python catmap-26.py '+reactant)
#print '='*20, 'removing C-C steps with TS > 1.5 eV'
#os.system('post-process-C-C-steps-v6.py '+reactant+' > catmap-rxns'+reactant)
#print 'did catmap remove --O assisted dehydrogenation'
print '='*20, 'graphviz'
os.system('python graphviz-17.py '+reactant)
#print '='*20, 'adsorption energy'
#os.system('python adsorption_energies-v42.py '+reactant)
#print '='*20, 'barriers'
#os.system('barriers-v13.py '+reactant)

