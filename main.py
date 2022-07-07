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

	# Decode the replacement format
	replacements = {
		'wild_type_AA': re.match(r'^[A-Z]', replacement),
		'position_of_AA': re.match(r'[0-9]+', replacement),
		'mutant_type_AA': re.match(r'[A-Z]$', replacement)
	}

	for element in replacements.values():
		element = element.group(0) if element is not None else element

	for types in AA_GROUPS.keys():
		for AA in AA_GROUPS[types]:
			pass
			# [TODO]
	
	return replacements

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
		if match:
			return False
		else:
			return True

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

	print("Groups filter: " + str(groups_check(INPUT, 'A345G')))
	print("Dipeptides?: " + str(dipeptide_matches(INPUT)))
	print("Inteins?: " + str(intein_matches(INPUT)))

#----------------------------------------
__main__()
