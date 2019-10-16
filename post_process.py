import os,sys
from copy import deepcopy
import pickle

def post_process(reactant):
	'open the file and make all the lines in one list'
	f = open(reactant+'-catmap.log')
	lines = f.read().split("\n")

	'a list of all elementary steps, including duplicates'
	elementary_steps = []
	for line in lines:
		elementary_steps.append(line.split('<->'))

	'getting information about the first elementary step to pre-initiate the loop'
	previous_step_reactant = str(elementary_steps[0]).split('+')[0]
	previous_step_reactant = previous_step_reactant[3:len(previous_step_reactant)-3]
	fs = []  #it should be empty before start. Otherwise, first step would be considered a duplicate
	post_process_catmap = []
	k,remove = 0,[] #list of numbers of steps to be removed

	'adding non-duplicate steps to a post-process list'
	for step in elementary_steps:
	  k += 1
	  if step == ['']:
  		'dealing with usually the last step where it is empty'
	  	continue
	  else:
		'define the reactant in the elemntary step'
		step_reactant = step[0].split('+')[0]
		step_reactant = step_reactant[1:len(step_reactant)-3]
		if step_reactant == previous_step_reactant:
			'if we are still with the same reactant, check if it is a duplicate'
			if step[2] in fs:
				#print 'duplicated step', step[2]
				remove.append(k)
			else:
				fs.append(step[2])
				post_process_catmap.append(step)
		else:
			'if not, we need to update reactant and create a new empty fs'
			fs = []
			previous_step_reactant = step_reactant
			fs.append(step[2])
			post_process_catmap.append(step)

	new = ''
	for i in remove:
		new = new+str(i)+'d;'
	os.system("sed '"+new+"' < "+reactant+"-catmap.log > "+reactant+"-catmap-post-process.log")
