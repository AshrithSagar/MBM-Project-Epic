"""
Peptide modifications.
"""
#========================================
import random
import argparse
import re
import csv
import pandas as pd
from scipy.stats import norm

#========================================
# DDG values from Alanine scan.
#----------------------------------------
def ddg_values_bals(alascan_file):
	"""Accepts .csv files: {Position, AA, DDG}."""
	df = pd.read_csv(alascan_file)
	df.set_index(['Index'])

	ddg_array = df[['Number', 'Name', 'IntraDDG']]
	# print("DDG values:\n", ddg_array)

	stable_ddg_values = ddg_array[ddg_array['IntraDDG'] < 0]
	print("Stable DDG values:\n", stable_ddg_values)

	return stable_ddg_values

def ddg_values_replot(alascan_file):
	"""Accepts .csv files from replot."""
	df = pd.read_csv(alascan_file)

	ddg_array = df[['ResNumber', 'ResName', 'ddGs']]
	# print("DDG values:\n", ddg_array)

	pos_ddg_values = ddg_array[ddg_array['ddGs'] > 0]
	print("Positive DDG values:\n", pos_ddg_values.sort_values('ddGs', ascending=False))

	neg_ddg_values = ddg_array[ddg_array['ddGs'] < 0]
	print("Negative DDG values:\n", neg_ddg_values.sort_values('ddGs', ascending=True))

	print("Ascending DDG values:\n", ddg_array.sort_values('ddGs', ascending=True))

	less_one_ddg_values = ddg_array[ddg_array['ddGs'] < 1]
	print("Less than threshold 1kcal/mol DDG values:\n", less_one_ddg_values.sort_values('ddGs', ascending=True))

	return less_one_ddg_values

def ddg_get_positions(dataframe):
	"""Accepts pandas dataframe"""
	df_ResNumber = dataframe[['ResNumber']]
	positions = []

	for index, ResNumber in df_ResNumber.iterrows():
		position = re.search(r'([0-9])+', str(ResNumber))
		position = position[0] if position is not None else position
		positions.append(position)
	return positions

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
	print("Original sequence:", sequence)
	mutated_sequences = []

	for mutation in mutations:
		print("Mutation:", mutation)
		# Recognise the group of the mutation.
		for types in AA_GROUPS.keys():
			if mutation['wild_type'] in AA_GROUPS[types]:
				print("The mutation is of type:", types)
				for AA in AA_GROUPS[types]:
					if not AA is mutation['wild_type']:
						# Replace the mutation with all possibilities within the group
						position = int(mutation['position'])
						print("Mutating", mutation['wild_type'],
							"with", AA, "at", position)
						new_sequence = sequence[:position-1] + AA + sequence[position:]
						mutated_sequences.append(new_sequence)
				break
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
def randomise_mutater(sequence):
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
# Random sampling through random, and/or Monte-Carlo.
#----------------------------------------
def random_sampler(sequences, approach, choose = 5):
	"""Random sampling"""
	if approach == 'Monte-Carlo':
		new_sequences = []
		# mean, var, skew, kurt = norm.stats(moments='mvsk')
		# [TODO]
		return new_sequences

	if approach == 'random':
		# Choose 5 sequences by default.
		new_sequences = random.sample(sequences, choose)
		return new_sequences

#========================================s

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
		ddg_preferences = ddg_values_replot(args.alaninescan)
		mutation_positions = ddg_get_positions(ddg_preferences)
		contents = [sequence]
		contents.extend(mutation_positions)
		print(contents)
		sequence, mutations = format_input(contents)
		sequences = groups_mutations(sequence, mutations)

	# sequences = intein_matches(sequences)

	print("Output sequences:", sequences)

	# If --output not specified, use input_file filename.
	output_file = args.output_file if args.output_file else args.input_file.replace(".txt", "_sequences.txt")
	
	# Save the output.
	new_contents = "\n".join(sequences)
	file = open(output_file, "w")
	file.write(new_contents)


if __name__ == "__main__":
	_main()
