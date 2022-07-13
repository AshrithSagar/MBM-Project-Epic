"""
Peptide modifications.
"""
#========================================
import random
import argparse
import re
from scipy.stats import norm

#========================================
# DDG values from Alanine scan.
#----------------------------------------
def ddg_values(alascan_file):
	"""Accepts .csv files: {Mutation, Position, DDG}."""
	with open(alascan_file, "r") as file:
	    contents = file.readlines()

	

#========================================
# Group mutations filter.
#----------------------------------------
def groups_mutations(sequence, mutations):
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
		# Recognise the group of the mutation.
		for types in AA_GROUPS.keys():
			if replacement['wild_type'] in AA_GROUPS[types]:
				print("The mutation is of type:", types)
				for AA in AA_GROUPS[types]:
					if not AA is replacement['wild_type']:
						# Replace the mutation with all possibilities within the group
						position = int(replacement['position'])
						print("Mutating", replacement['wild_type'],
							"with", AA, "at", position)
						new_sequence = sequence[:position] + AA + sequence[position+1:]
						mutated_sequences.append(new_sequence)

		# Random sampling
		# Through Monte-Carlo.
		# mean, var, skew, kurt = norm.stats(moments='mvsk')
		

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
	# 	sequences = groups_mutations(sequence, [str(mutation_position)])
	
	# mutation_positions = range(0, len(sequence))
	# new_mutation_positions = []
	# for position in map(str, mutated_positions):
	# 	new_mutation_positions.append(str(position))
	# mutated_positions = new_mutation_positions
	# sequences = groups_mutations(sequence, mutation_positions)

#========================================
def format_input(contents):
	sequence = contents[0].replace('\n', '')
	given_mutations = contents[1:]
	mutations = []

	# Decode the replacement mutations format.
	for mutation in given_mutations:
		mutation = mutation.replace('\n', '') # Remove newlines.
		
		replacement = {
			'wild_type': re.search(r'^([A-Za-z])', mutation),
			'position': re.search(r'([0-9])+', mutation),
			'mutant_type': re.search(r'([A-Za-z])$', mutation) # Ignored currently.
		}

		for key in replacement.keys():
			element = replacement[key]
			replacement[key] = element[0] if element is not None else element

		# If only position is given.
		if replacement['wild_type'] == None:
			position = int(replacement['position'])
			replacement['wild_type'] = sequence[position]

		mutations.append(replacement)
	return sequence, mutations

def _main():
	# For the command line parser.
	parser = argparse.ArgumentParser(prog='pep_mod', description='Peptide modifications generator')
	parser.add_argument('input_file', type=str, help='Input file [.txt]')
	parser.add_argument('-o', '--output', dest='output_file', type=str, help='Output filename')
	parser.add_argument('-d', '--dipeptide', action='store_true', help='Dipeptides match')
	parser.add_argument('-g', '--groups', action='store_true', help='Groups filter')
	parser.add_argument('-a', '--alaninescan', help='Alanine scan DDG results')
	args = parser.parse_args()

	with open(args.input_file, "r") as file:
	    contents = file.readlines()
	sequence, mutations = format_input(contents)
	sequences = [sequence]

	if args.groups:
		sequences = groups_mutations(sequence, mutations)

	if args.dipeptide:
		sequences = dipeptide_matches(sequences)
	
	if args.alaninescan:
		sequences = ddg_values(args.alaninescan)

	# sequences = intein_matches(sequences)

	print("Output sequences:", sequences)

	# If --output not specified, use input_file filename.
	output_file = args.output_file if args.output_file else args.input_file.replace(".txt", "_sequences.txt")
	
	# Save the output.
	new_contents = "".join(sequences)
	file = open(output_file, "w")
	file.write(new_contents)


if __name__ == "__main__":
	_main()
