"""
Include pdb chains.
"""
import sys
import getopt

with open(sys.argv[1], "r") as file:
    contents = file.readlines()

chains = sys.argv[3]
new_file = ""

chain = chains.pop(0)
for line in contents:
	if line.startswith("TER"):
		chain = chains.pop(0)
		continue
	if line.startswith("ATOM"):
		new_line = line[:21] + new_chain + line[23:] + "\n"
	else:
		new_line = line
	try: new_file.append(new_line)
	except: pass
	
# Save as .json
file = sys.argv[2]
file = open(file, "w")
file.write(str(new_file))
