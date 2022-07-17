"""Peptide modifications
"""
import random
import argparse
import re
import csv
import itertools
from pandas import read_csv
from scipy.stats import norm


class BAlaS:
	"""BUDE Alanine Scan: ddG values.
	"""
	def __init__(self, alascan_file):
		self.stable_df_bals = self.df_bals = self.df_ddg = None
		self.positions = []

	def bals_read(self, bals2csv_file):
		"""Accepts .csv files converted using bals2csv.py
		"""
		self.df_bals = read_csv(bals2csv_file)
		self.df_bals.set_index(['Index'])

		self.df_bals = self.df_bals[['Number', 'Name', 'IntraDDG']]
		# print("DDG values:\n", self.ddg_bals)

		self.stable_df_bals = self.df_bals[self.df_bals['IntraDDG'] < 0]
		print("Stable DDG values:\n", self.stable_df_bals)


	def replot_read(self, replot_csv_file):
		"""Accepts .csv files from replot"""
		self.df_ddg = read_csv(replot_csv_file)
		self.df_ddg = self.df_ddg[['ResNumber', 'ResName', 'ddGs']]
		# print("DDG values:\n", self.df_ddg)


	def replot_filter(self, threshold):
		"""Filter out based on ddG thresholds"""
		df = self.df_ddg
		# ddg_pos = df[df['ddGs'] > 0]
		# print("Positive DDG values:\n", ddg_pos.sort_values('ddGs', ascending=False))

		# ddg_neg = df[df['ddGs'] < 0]
		# print("Negative DDG values:\n", ddg_neg.sort_values('ddGs', ascending=True))

		# print("Ascending DDG values:\n", df.sort_values('ddGs', ascending=True))

		ddg_threshold = df[df['ddGs'] < threshold]
		print("DDG values for less than threshold", threshold, ":\n",
			ddg_threshold.sort_values('ddGs', ascending=True))


	def replot_get_positions(self):
		"""Get positions from filtered ddG thresholds"""
		df_ResNumber = self.df_ddg[['ResNumber']]
		self.positions = []

		for index, ResNumber in df_ResNumber.iterrows():
			position = re.search(r'([0-9])+', str(ResNumber))
			position = position[0] if position is not None else position
			self.positions.append(position)


#========================================
# Group mutations filter.
#----------------------------------------
class mutations:
	def __init__(self):
		self.AA_GROUPS = {
			'polar_uncharged': ['S', 'T', 'C', 'N', 'Q'],
			'positively_charged': ['K', 'R', 'H'],
			'negatively_charged': ['D', 'E'],
			'nonpolar_aliphatic': ['G', 'A', 'V', 'L', 'M', 'I'],
			'nonpolar_aromatic': ['F', 'Y', 'W']
		}
		self.sequence = self.mutations = None


	def by_groups(self):
		"""Returns mutation arrays based on group mutations
		"""
		print("Original sequence:", sequence)
		new_mutations = []

		for mutation in mutations:
			print("Mutation:", mutation)
			if mutation in mutation_lock:
				print("G: Skipping locked position", mutation['position'])
				continue
			if mutation['wild_type'] == 'P':
				print("Skipping Proline mutation at", mutation['position'])
				continue
			if mutation['mutant_type']:
				AA = mutation['mutant_type']
				position = int(mutation['position'])
				print("Mutating", mutation['wild_type'],
					"with", AA, "at", position)
				new_mutations.append(mutation)
			else:
				# Recognise the group of the mutation.
				for types in AA_GROUPS.keys():
					if mutation['wild_type'] in AA_GROUPS[types]:
						# Self mutations are included.
						print("The mutation is of type:", types)
						for AA in AA_GROUPS[types]:
							# Replace the mutation with all possibilities within the group
							position = int(mutation['position'])
							print("G: Mutating", mutation['wild_type'],
								"with", AA, "at", position)
							new_mutation = mutation['wild_type'] + mutation['position'] + AA
							new_mutations.append(new_mutation)
						break
		return new_mutations

	def to_sequences(self):
		"""P&C of new_mutations. nCr approach.
		Choose r mutation positions at a time, out of n mutations"""
		# sequence, mutations, count
		permutations = itertools.permutations(mutations, count)

		mutated_sequences = []
		for mutation_tuple in permutations:
			print("S: Mutations:", mutation_tuple)
			for mutation in mutation_tuple:
				replacement = {
					'wild_type': re.search(r'^([A-Za-z])', mutation),
					'position': re.search(r'([0-9])+', mutation),
					'mutant_type': re.search(r'([A-Za-z])$', mutation)
				}

				for key in replacement.keys():
					element = replacement[key]
					replacement[key] = element[0] if element is not None else element

				position = int(replacement['position'])
				new_sequence = sequence[:position-1] + replacement['mutant_type'] + sequence[position:]
				mutated_sequences.append(new_sequence)

		return mutated_sequences

	def dipeptide_match(sequence):
		"""Returns whether a dipeptide is present or not"""
		match = re.search(r"(.)\1", str(sequence)) # Returns first dipeptide match.
		if match: return match[0], match.span()
		return None, None

	def dipeptide_matches(sequences):
		"""Filters out and returns the non-dipeptide sequences"""
		new_sequences = sequences
		for sequence in sequences:
			if dipeptide_match(sequence):
				new_sequences.pop(sequence)

		print("Dipeptides filtered: " + str(new_sequences))
		return new_sequences

	def dipeptide_mutater(sequence, dipeptide, ddg):
		"""Mutates to remove dipeptides."""
		AA_GROUPS = {
			'polar_uncharged': ['S', 'T', 'C', 'P', 'N', 'Q'],
			'positively_charged': ['K', 'R', 'H'],
			'negatively_charged': ['D', 'E'],
			'nonpolar_aliphatic': ['G', 'A', 'V', 'L', 'M', 'J', 'I'],
			'nonpolar_aromatic': ['F', 'Y', 'W']
		}
		match, span = dipeptide[0], dipeptide.span()

		# Decide mutation position by comparing ddG values.
		position = span[0]
		lesser_ddg = (ddg['ddGs'][position] > ddg['ddGs'][position+1])
		mutation_position = position+2 if lesser_ddg else position+1
		print("Position", mutation_position, "has lesser ddG value among the dipeptide", match)

		# Conservative replacement.
		contents = [sequence]
		contents.append(str(mutation_position)) # Just one mutation_position.
		print(contents)
		sequence, mutations = format_input(contents)

		# Discard mutations that produce another dipeptide.
		new_mutations = []
		for mutation in mutations:
			if mutation in mutation_lock:
				print("DP: Skipping locked position", mutation['position'])
				continue
			# Recognise the group of the mutation.
			for types in AA_GROUPS.keys():
				if mutation['wild_type'] in AA_GROUPS[types]:
					for AA in AA_GROUPS[types]:
						position = int(mutation['position'])
						if (AA == sequence[position-2])|(AA == sequence[position]):
							print("Discarding mutation of", mutation['wild_type'],
								"with", AA, "at", position)
						else:
							print("DP: Mutating", mutation['wild_type'],
								"with", AA, "at", position)
							new_mutation = mutation['wild_type'] + mutation['position'] + AA
							new_mutations.append(new_mutation)
					break
		return new_mutations
	

	def by_intein_sequences(self):
		"""[SKIPPED]
		"""
		pass


	def by_cleavage_sites(self):
		"""[SKIPPED]
		"""
		pass
		# def intein_matches(sequence):
		# inteins = ["CRAZY_SEQUENCE", "ANOTHER_SEQUENCE"]
		# # Put the source as a dictionary or an array, preferably.

		# print("Inteins?: " + str(intein_matches(sequence)))
		# for intein in inteins:
		# 	match = str(sequence).find(intein)
		# 	return not match


	def by_charge_criterion(self):
		"""[SKIPPED]
		"""
		pass


	def randomise(self):
		"""At random positions
		"""
		pass
		# [TODO]
		# for mutation_position in range(len(sequence)):
		# 	sequences = groups_mutations(sequence, [str(mutation_position)])

		# mutation_positions = range(0, len(sequence))
		# new_mutation_positions = []
		# for position in map(str, mutated_positions):
		# 	new_mutation_positions.append(str(position))
		# mutated_positions = new_mutation_positions
		# sequences = groups_mutations(sequence, mutation_positions)


	def random_sampler(self, choose=5):
		"""Random sampling through random.
		Chooses 5 sequences by default.
		"""
		new_sequences = random.sample(sequences, choose)
		return new_sequences


	def monte_carlo_sampler(self, choose=5):
		"""Random sampling through Monte-Carlo.
		Chooses 5 sequences by default.
		"""
		new_sequences = []
		# mean, var, skew, kurt = norm.stats(moments='mvsk')
		# [TODO]
		return new_sequences


def save_sequences(file, sequences):
	contents = "\n".join(sequences)
	file = open(file, "w")
	file.write(contents)


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
			'mutant_type': re.search(r'([A-Za-z])$', mutation)
		}

		for key in replacement.keys():
			element = replacement[key]
			replacement[key] = element[0] if element is not None else element

		# If only position is given.
		if replacement['wild_type'] == None:
			position = int(replacement['position'])
			replacement['wild_type'] = sequence[position-1]

		mutations.append(replacement)
	return sequence, mutations


def main():
	# For the command line parser.
	parser = argparse.ArgumentParser(prog='pep_mod', description='Peptide modifications generator')
	parser.add_argument('input_file', type=str, help='Input file [.txt]')
	parser.add_argument('-o', '--output', dest='output_file', type=str, help='Output filename')
	parser.add_argument('-d', '--dipeptide', action='store', help='Dipeptides match')
	parser.add_argument('-g', '--groups', action='store_true', help='Groups filter')
	parser.add_argument('-a', '--alaninescan', help='Alanine scan DDG results')
	parser.add_argument('-l', '--lock', type=str, dest='mutation_lock', help='Mutation lock positions')
	parser.add_argument('-c', '--count', type=int, dest='mutation_count', help='Number of mutations to consider at a time', default=1)
	args = parser.parse_args()

	with open(args.input_file, "r") as file:
	    input_contents = file.readlines()
	sequence = input_contents[0]
	sequences = [sequence]

	if args.mutation_lock:
		with open(args.mutation_lock, "r") as file:
			lock_contents = file.readlines()
		global mutation_lock
		mutation_lock = [sequence]
		for line in lock_contents:
			content = line.replace('\n', '') # Remove newlines.
			mutation_lock.append(content)
		print("Locked mutation positions:", mutation_lock[1:])
		sequence, mutation_lock = format_input(mutation_lock)

	if args.groups:
		sequence, mutations = format_input(input_contents)
		sequences = groups_mutations(sequence, mutations)

	if args.alaninescan:
		ddg_array = ddg_replot_read(args.alaninescan)
		ddg_preferences = ddg_replot_filter(ddg_array)
		mutation_positions = ddg_get_positions(ddg_preferences)
		contents = [sequence]
		contents.extend(mutation_positions)
		print(contents)
		sequence, mutations = format_input(contents)
		mutations = groups_mutations(sequence, mutations)
		sequences = mutations2sequences(sequence, mutations, args.mutation_count)
		print("="*50)

	if args.dipeptide:
		# sequences = dipeptide_matches(sequence)
		for seq in sequences:
			contents = [seq]
			ddg_array = ddg_replot_read(args.dipeptide)
			matches = re.finditer(r"(.)\1", str(seq)) # Find all dipeptides.
			for match in matches:
				print(match)
				mutations = dipeptide_mutater(seq, match, ddg_array)
				contents.extend(mutations)
			print(contents)
			seq, all_mutations = format_input(contents)
			seqs = groups_mutations(seq, all_mutations)
		sequences.extend(seqs)

	# sequences = intein_matches(seqs)

	# Get unique sequences.
	sequences = list(dict.fromkeys(sequences))

	print("Output sequences:", sequences)
	print("Output sequences count:", len(sequences))

	# If --output not specified, use input_file filename.
	output_file = args.output_file if args.output_file else args.input_file.replace(".txt", "_allSequences.txt")
	save_sequences(output_file, sequences)

	try:
		sequences = random_sampler(sequences, 'random', 5)
		output_file = args.output_file if args.output_file else args.input_file.replace(".txt", "_randomSequences.txt")
		save_sequences(output_file, sequences)
	except: pass


if __name__ == "__main__":
	main()
