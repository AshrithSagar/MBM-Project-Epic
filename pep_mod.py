"""Peptide modifications
"""
import random
import argparse
import re
import csv
import copy
import itertools
import pandas
from scipy.stats import norm


class BAlaS:
	"""BUDE Alanine Scan: ddG values
	"""
	def __init__(self, alascan_file):
		self.stable_df_bals = self.df_bals = self.df_ddg = None
		self.positions = []

	def bals_read(self, bals2csv_file):
		"""Accepts .csv files converted using bals2csv.py
		"""
		self.df_bals = pandas.read_csv(bals2csv_file)
		self.df_bals.set_index(['Index'])

		self.df_bals = self.df_bals[['Number', 'Name', 'IntraDDG']]
		# print("DDG values:\n", self.ddg_bals)

		self.stable_df_bals = self.df_bals[self.df_bals['IntraDDG'] < 0]
		print("Stable DDG values:\n", self.stable_df_bals)


	def replot_read(self, replot_csv_file):
		"""Accepts .csv files from replot"""
		self.df_ddg = pandas.read_csv(replot_csv_file)
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
			if position is not None: position = int(position[0])
			self.positions.append(position)


class mutation_object:
	"""Mutation format
	"""
	def __init__(self, sequence, mutation):
		self.wild_type = self.position = self.mutant_type = None
		self.sequence = sequence
		self.from_str(mutation)


	def from_str(self, mutation):
		position = re.search(r'([0-9])+', mutation)
		self.position = int(position[0]) if position is not None else self.position

		wild_type = re.search(r'^([A-Za-z])', mutation)
		self.wild_type = wild_type[0] if wild_type is not None else self.wild_type
		if self.wild_type is None: self.wild_type = self.sequence[self.position-1]

		mutant_type = re.search(r'([A-Za-z])$', mutation)
		self.mutant_type = mutant_type[0] if mutant_type is not None else self.mutant_type


	def new_mutant_type(self, mut_ty):
		new_obj = copy.copy(self)
		new_obj.mutant_type = mut_ty
		print("->", vars(self))
		return new_obj


	def to_str(self):
		mutant_type = self.mutant_type if self.mutant_type else ''
		mut_str = self.wild_type + str(self.position) + mutant_type
		return mut_str


class mutater:
	"""Mutater: groups, ddg, dipeptide
	"""
	def __init__(self, sequence, mutations=[], mutation_lock=[]):
		self.AA_GROUPS = {
			'polar_uncharged': ['S', 'T', 'C', 'N', 'Q'],
			'positively_charged': ['K', 'R', 'H'],
			'negatively_charged': ['D', 'E'],
			'nonpolar_aliphatic': ['G', 'A', 'V', 'L', 'M', 'I'],
			'nonpolar_aromatic': ['F', 'Y', 'W']
		}
		self.sequence = sequence
		self.mutations = mutations
		self.mutation_lock = mutation_lock


	def append_mutation(self, mut_obj):
		"""Append input mutation object to self.mutations
		"""
		self.mutations.append(mut_obj)


	def by_groups(self):
		"""Returns mutation arrays based on group mutations
		"""
		sequence = self.sequence
		mutation_lock = self.mutation_lock
		groups = self.AA_GROUPS
		print("G: Original sequence:", sequence)
			
		new_mutations = []
		for mutation in self.mutations:
			print("G: Mutation:", mutation.to_str())
			if mutation in mutation_lock:
				print("G: Skipping locked position", mutation.position)
				continue
			if mutation.wild_type == 'P':
				print("G: Skipping Proline mutation at", mutation.position)
				continue
			if mutation.mutant_type:
				print("G: Mutating", mutation.wild_type,
					"with", mutation.mutant_type, "at", mutation.position)
				new_mutations.append(mutation)
			else:
				# Recognise the group of the mutation.
				for types in groups.keys():
					if mutation.wild_type in groups[types]:
						# Self mutations are included.
						print("G: The mutation is of type:", types)
						for AA in groups[types]:
							# Replace the mutation with all possibilities within the group
							print("G: Mutating", mutation.wild_type,
								"with", AA, "at", mutation.position)
							new_mutation = mutation.new_mutant_type(AA)
							new_mutations.append(new_mutation)
						break
		self.mutations = new_mutations
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
		A normal distribution is chosen.
		Chooses 5 sequences by default.
		"""
		new_sequences = []
		# mean, var, skew, kurt = norm.stats(moments='mvsk')
		# [TODO]
		return new_sequences


class sequences:
	def __init__(self):
		pass

	def check_dipeptide(self, sequence):
		"""Returns first dipeptide match if present, else None
		"""
		match = re.search(r"(.)\1", str(sequence))
		if match: return match[0], match.span()
		return None, None

	def dipeptide_matches(self, sequences):
		"""Filters out and returns the non-dipeptide sequences"""
		new_sequences = sequences
		for sequence in sequences:
			if dipeptide_match(sequence):
				new_sequences.pop(sequence)

		print("Dipeptides filtered: " + str(new_sequences))
		return new_sequences

	def dipeptide_mutater(self, sequence, dipeptide, ddg):
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


def save_as(file, contents):
	content = "\n".join(contents)
	file = open(file, "w")
	file.write(content)


def format_input(contents):
	remove_new_line = lambda s: s.replace('\n', '')
	
	sequence = remove_new_line(contents[0])
	given_mutations = contents[1:]
	mutations_obj = mutater(sequence)
	
	# Decode in the mutations format.
	for line in given_mutations:
		line = remove_new_line(line)
		mut = mutation_object(sequence, line)
		mutations_obj.append_mutation(mut)
	return mutations_obj


def main():
	# For the command line parser.
	parser = argparse.ArgumentParser(prog='pep_mod', description='Peptide modifications generator')
	parser.add_argument('input_file', type=str, help='Input file [.txt]')
	parser.add_argument('-o', '--output', dest='output_file', type=str, help='Output filename')
	parser.add_argument('-d', '--dipeptide', action='store', help='Dipeptides match')
	parser.add_argument('-g', '--groups', action='store_true', help='Groups filter')
	parser.add_argument('-a', '--alaninescan', help='Alanine scan DDG results from BUDE Alanine scan')
	parser.add_argument('-l', '--lock', type=str, dest='mutation_lock', help='Mutation lock positions')
	parser.add_argument('-c', '--count', type=int, dest='mutation_count', help='Number of mutations to consider at a time', default=1)
	args = parser.parse_args()

	# If --output not specified, use input_file filename.
	output_file = args.output_file if args.output_file else args.input_file

	with open(args.input_file, "r") as file:
	    input_contents = file.readlines()
	mutations_obj = format_input(input_contents)
	sequence = mutations_obj.sequence
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
		mutations = mutations_obj.by_groups()
		for mutation in mutations:
			print(mutation.to_str())
		exit(0)
		sequences = mutate.to_sequences(mutations, args.mutation_count)
		sequence, mutations = format_input(input_contents)
		sequences = groups_mutations(sequence, mutations)

	if args.alaninescan:
		bude = BAlaS()
		df_ddg = bude.replot_read(args.alaninescan)
		ddg_preferences = bude.replot_filter(df_ddg, 1) # Threshold 1kCal/mol
		positions = bude.get_positions(ddg_preferences)
		# [TODO]
		mutate = mutations(sequence)
		mutations = mutate.by_groups()
		sequences = mutate.to_sequences(mutations, args.mutation_count)
		print("="*50)

	if args.dipeptide:
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

	get_unique = lambda seqs: list(dict.fromkeys(seqs))
	sequences = get_unique(sequences)

	print("Output sequences:", sequences)
	print("Output sequences count:", len(sequences))

	save_as(output_file.replace(".txt", "_allSequences.txt"), sequences)

	try:
		sequences = random_sampler(sequences, 'random', 5)
		save_as(output_file.replace(".txt", "_randomSequences.txt"), sequences)
	except: pass


if __name__ == "__main__":
	main()
