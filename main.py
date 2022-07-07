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
			'wild_type_AA': re.findall(r'^(\w)', mutation)[0],
			'position_of_AA': re.findall(r'(\d)+', mutation)[0],
			'mutant_type_AA': re.findall(r'(\w)$', mutation)[0] # Ignored currently.
		}

		# Replace the mutation with all possibilities within the group,
		# by recognising its group.
		for types in AA_GROUPS.keys():
			if replacement['wild_type_AA'] in AA_GROUPS[types]:
				print("The mutation is of type:", types)
				for AA in AA_GROUPS[types]:
					if not AA is replacement['wild_type_AA']:
						print("Mutating", replacement['wild_type_AA'], "with", AA)
						new_sequence = sequence
						new_sequence.replace(AA, replacement['position_of_AA'])
						mutated_sequences.append(new_sequence)

		# Random sampling
			# Through Monte-Carlo.
				# [TODO]

			# The approach using random.
		# mutated_sequences = random.shuffle(mutated_sequences)[0:5]

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
	return new_sequences

#========================================
# Intein sequences.
#----------------------------------------
# SKIPPED, FOR NOW.
def intein_matches(sequence):
	inteins = ["CRAZY_SEQUENCE", "ANOTHER_SEQUENCE"] 
	# Put the source as a dictionary or an array, preferably.

	for intein in inteins:
		match = str(sequence).find(intein)
		return not match

#========================================
# Cleavage sites.
#----------------------------------------
# [TODO] Possibly skipped, for now.

#========================================
# Charge criterion
#----------------------------------------
# [TODO]

#========================================
def __main__():
	INPUT = "ERCYDNTAGTSYVVGETWEKPYQGWMIVDCTCLGEGSGRITCT" # Peptide sequence here.

	print("Groups filtered: " + str(groups_check(INPUT, 
		['R2', 'W20'] # Enter the mutations required here.
		)))
	print("Dipeptides filtered: " + str(dipeptide_matches(INPUT)))
	print("Inteins?: " + str(intein_matches(INPUT)))

#----------------------------------------
__main__()
