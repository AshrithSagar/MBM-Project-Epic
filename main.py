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
def groups_check(sequence, replacement):
	AA_GROUPS = {
		'polar_uncharged': ['S', 'T', 'C', 'P', 'N', 'Q'],
		'positively_charged': ['K', 'R', 'H'],
		'negatively_charged': ['D', 'E'],
		'nonpolar_aliphatic': ['G', 'A', 'V', 'L', 'M', 'J', 'I'],
		'nonpolar_aromatic': ['F', 'Y', 'W']
	}
	sequence = str(sequence)
	replacement = str(replacement)
	sequences = [] # The input for the next iteration.

	# Decode the replacement format.
	replacements = {
		'wild_type_AA': re.match(r'^[A-Z]', replacement),
		'position_of_AA': re.match(r'[0-9]+', replacement),
		'mutant_type_AA': re.match(r'[A-Z]$', replacement)
	}
	for element in replacements.values():
		element = element.group(0) if element is not None else element

	# Replace the mutation with all possibilities within the group, by recognising its group.
	for types in AA_GROUPS.keys():
		if replacements['mutant_type_AA'] in AA_GROUPS[types]:
			for AA in AA_GROUPS[types]:
				print(AA)
				new_sequence = sequence.copy()
				new_sequence[replacements['position_of_AA']] = AA
				sequences.append(new_sequence)

	# Random sampling through Monte-Carlo.
	# [TODO]

	return sequences

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
	inteins = ["CRAZY_SEQUENCE", "ANOTHER_SEQUENCE"] # Put the source as a dictionary or an array, preferably.

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

	print("Groups filter: " + str(groups_check(INPUT, 'A2G')))
	print("Dipeptides?: " + str(dipeptide_matches(INPUT)))
	print("Inteins?: " + str(intein_matches(INPUT)))

#----------------------------------------
__main__()
