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
# AA_GROUPS = {
# 	'aliphatic' : ["G", "A", "V", "L", "I", "P"]
# 	'polar neutral' : ["S", "T"],
# 	'amide' : ["N", "Q"],
# 	'sulfur-containing' : ["C", "M"],
# 	'aromatic' : ["F", "Y", "W"],
# 	'anionic' : ["D", "E"],
# 	'cationic' : ["H", "K", "R"],
# 	'hydrophobic' : ["V", "I", "L", "F", "W", "Y", "M"],
# 	'hydrophilic' : ["S", "T", "H", "N", "Q", "E", "D", "K", "R"],
# 	'positively-charged' : ["K", "R", "H"],
# 	'negatively-charged' : ["D", "E"]
# }

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
	SOURCE_GROUP = {} # Put the source as a dictionary, preferably.

#========================================
# Cleavage sites.
#----------------------------------------
# [TODO]

#========================================
# Charge criterion
#----------------------------------------
# [TODO]

#========================================
def __main__():
	INPUT = "HHHDAASDLKJASD" # Peptide sequence here.
	
	print(dipeptide_matches(INPUT))

#----------------------------------------
__main__()
