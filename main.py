"""
Peptide modifications.
"""
#========================================
import random
import re
import numpy
import scipy

#========================================
# Alanine scan.
#----------------------------------------
# [TODO]

#========================================
# Groups
def groups_check(sequence, mutations):
	AA_GROUPS = {
		'polar_uncharged': ['S', 'T', 'C', 'P', 'N', 'Q'],
		'positively_charged': ['K', 'R', 'H'],
		'negatively_charged': ['D', 'E'],
		'nonpolar_aliphatic': ['G', 'A', 'V', 'L', 'M', 'J', 'I'],
		'nonpolar_aromatic': ['F', 'Y', 'W']
	}
	sequence = str(sequence)
	print("Original sequence:", sequence)
	mutated_sequences = [] # The inputs for the next iteration.

	for mutation in mutations:
		mutation = str(mutation)
		print("Mutation:", mutation)

		# Decode the replacement format.
		replacement = {
			'wild_type': re.findall(r'^(\w)', mutation)[0],
			'position': re.findall(r'(\d)+', mutation)[0],
			'mutant_type': re.findall(r'(\w)$', mutation)[0] # Ignored currently.
		}
		if replacement['wild_type'] == None:
			position = int(replacement['position'])
			replacement['wild_type'] = sequence[position]

		# Replace the mutation with all possibilities within the group,
		# by recognising its group.
		for types in AA_GROUPS.keys():
			if replacement['wild_type'] in AA_GROUPS[types]:
				print("The mutation is of type:", types)
				for AA in AA_GROUPS[types]:
					if not AA is replacement['wild_type']:
						position = int(replacement['position'])
						print("Mutating", replacement['wild_type'],
							"with", AA, "at", position)
						new_sequence = sequence[:position] + AA + sequence[position+1:]
						mutated_sequences.append(new_sequence)

		# Random sampling
		# Through Monte-Carlo.
		# [TODO]

		# The approach using random.
		mutated_sequences = random.sample(mutated_sequences, 5)

	print("Groups filtered: " + str(mutated_sequences))
	return mutated_sequences

#========================================
# Dipeptide filter.
#----------------------------------------
def dipeptide_match(sequence):
	"""Returns whether a dipeptide is present or not"""
	matches = re.search(r"(.)\1", str(sequence))
	return bool(matches)

def dipeptide_matches(sequences):
	"""Filters out and returns the non-dipeptide sequences"""
	new_sequences = sequences
	for sequence in sequences:
		if dipeptide_match(sequence):
			new_sequences.pop(sequence)

	print("Dipeptides filtered: " + str(new_sequences))
	return new_sequences

#========================================
# Intein sequences.
#----------------------------------------
# SKIPPED, FOR NOW.
# def intein_matches(sequence):
# 	inteins = ["CRAZY_SEQUENCE", "ANOTHER_SEQUENCE"] 
# 	# Put the source as a dictionary or an array, preferably.

# 	print("Inteins?: " + str(intein_matches(sequence)))
# 	for intein in inteins:
# 		match = str(sequence).find(intein)
# 		return not match

#========================================
# Cleavage sites.
#----------------------------------------
# [TODO] Possibly skipped, for now.

#========================================
# Charge criterion
#----------------------------------------
# [TODO]

#========================================
# Randomiser for positions.
#----------------------------------------
def randomise_position(sequence):
	"""Random positions"""
	pass
	# for mutation_position in range(len(sequence)):
	# 	sequences = groups_check(sequence, [str(mutation_position)])
	
	# mutation_positions = range(0, len(sequence))
	# new_mutation_positions = []
	# for position in map(str, mutated_positions):
	# 	new_mutation_positions.append(str(position))
	# mutated_positions = new_mutation_positions
	# sequences = groups_check(sequence, mutation_positions)

#========================================
def __main__():
	sequence = "ERCYDNTAGTSYVVGETWEKPYQGWMIVDCTCLGEGSGRITCT" # Peptide sequence here.

	sequences = groups_check(sequence, 
		['R2', 'W20'] # Enter the mutations required here.
	)
	# sequences = dipeptide_matches(sequences)
	# sequences = intein_matches(sequences)
	print("Output sequences:", sequences)

#----------------------------------------
__main__()
