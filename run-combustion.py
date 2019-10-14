import os
import sys

reactant = sys.argv[1]	#name of the hydrocarbon reactant

print 'generating reactions list ...'
os.system('python generate-species.py '+reactant)
print 'eliminating repeated reactions ...'
os.system('python reduce-species.py '+reactant) 
print 'writing reactions in catmap format ...'
os.system('python catmap.py '+reactant)
#print '='*20, 'graphviz'
#os.system('python graphviz-17.py '+reactant)

