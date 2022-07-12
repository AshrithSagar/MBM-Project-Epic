"""
Include pdb chains.
"""
import sys
import argparse

def _main():
	parser = argparse.ArgumentParser(description='Modify pdb chains')
	parser.add_argument('input_file', type=str)
	parser.add_argument('chains', type=list)
	parser.add_argument('-o', '--output', dest='outputfile', type=str, help='output filename')
	args = parser.parse_args()

	with open(args.input_file, "r") as file:
	    contents = file.readlines()

	print (args.chains)
	# new_file = ""

	# chain = chains.pop(0)
	# for line in contents:
	# 	if line.startswith("TER"):
	# 		chain = chains.pop(0)
	# 		continue
	# 	if line.startswith("ATOM"):
	# 		new_line = line[:21] + new_chain + line[23:] + "\n"
	# 	else:
	# 		new_line = line
	# 	try: new_file.append(new_line)
	# 	except: pass
		
	# # Save as .json
	# file = argv[2]
	# file = open(file, "w")
	# file.write(str(new_file))

if __name__ == "__main__":
	_main()
