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
# 

#========================================
# Groups
AA_GROUPS = {
	'aliphatic' : ["G", "A", "V", "L", "I", "P"],
	'polar neutral' : ["S", "T"],
	'amide' : ["N", "Q"],
	'sulfur-containing' : ["C", "M"],
	'aromatic' : ["F", "Y", "W"],
	'anionic' : ["D", "E"],
	'cationic' : ["H", "K", "R"],
	'hydrophobic' : ["V", "I", "L", "F", "W", "Y", "M"],
	'hydrophilic' : ["S", "T", "H", "N", "Q", "E", "D", "K", "R"],
	'positively-charged' : ["K", "R", "H"],
	'negatively-charged' : ["D", "E"]
}

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

	print("Dipeptides?: " + str(dipeptide_matches(INPUT)))
	print("Inteins?: " + str(intein_matches(INPUT)))

#----------------------------------------
__main__()
