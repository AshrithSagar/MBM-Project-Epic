"""
Include pdb chains.
"""
import sys
import getopt

def help():
	print("Usage: ", argv[0], '-i <inputfile> -o <outputfile>')

def main(argv):
	inputfile = ''
	outputfile = ''
	try:
		opts, args = getopt.getopt(argv,"hi:o:",["input=","output="])
		print(opts)
		print(args)
	except getopt.GetoptError:
		help()
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			help()
			sys.exit()
		elif opt in ("-i", "--input"):
			inputfile = arg
		elif opt in ("-o", "--output"):
		 	outputfile = arg
	print('Input file is ', inputfile)
	print('Output file is ', outputfile)

	# with open(argv[1], "r") as file:
	#     contents = file.readlines()

	# chains = argv[3]
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
	main(sys.argv)
