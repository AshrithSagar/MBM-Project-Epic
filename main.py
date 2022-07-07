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
	mutated_sequences = [] # The inputs for the next iteration.

	for mutation in mutations:
		mutation = str(mutation)

		# Decode the replacement format.
		replacement = {
			'wild_type_AA': re.match(r'^[A-Z]', mutation),
			'position_of_AA': re.match(r'[0-9]+', mutation),
			'mutant_type_AA': re.match(r'[A-Z]$', mutation)
		}
		for element in replacement.values():
			element = element.group(0) if element is not None else element
		print(replacement)

		# Replace the mutation with all possibilities within the group,
		# by recognising its group.
		for types in AA_GROUPS.keys():
			if replacement['mutant_type_AA'] in AA_GROUPS[types]:
				print("The mutation is of type: ", types)
				for AA in AA_GROUPS[types]:
					print("Mutating ", replacement['wild_type_AA'], "with", AA)
					new_sequence = sequence.copy()
					new_sequence[replacement['position_of_AA']] = AA
					mutated_sequences.append(new_sequence)

		# Random sampling through Monte-Carlo.
		# [TODO]

	return mutated_sequences

#========================================
# Dipeptide filter.
#----------------------------------------
def dipeptide_matches(sequence):
	"""Returns whether a dipeptide is present or not"""
	matches = re.search(r"(.)\1", str(sequence))
	return matches

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
	INPUT = "HHHDAASDLKJASD" # Peptide sequence here.

	print("Groups filter: " + str(groups_check(INPUT, 
		['R2K', 'A7G'] # Enter the mutations required here.
		)))
	print("Dipeptides?: " + str(dipeptide_matches(INPUT)))
	print("Inteins?: " + str(intein_matches(INPUT)))

#----------------------------------------
__main__()
