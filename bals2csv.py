"""
.bals to .csv converter
"""
import sys
import json
import pandas as pd

with open(sys.argv[1], "r") as file:
    contents = file.readlines()

json_file = {
	"Index": [],
	"Number": [],
	"Name": [],
	"Chain": [],
	"InterDG": [],
	"InterDDG": [],
	"NormTerDDG": [],
	"IntraDG": [],
	"IntraDDG": [],
	"NormTraDDG": [],
	"ChainAtoms": []
}

for line in contents:
	if line.startswith("#"): continue

	value = line[0:6].strip()
	try: json_file["Index"].append(int(value))
	except: pass

	value = line[6:13].strip()
	try: json_file["Number"].append(int(value))
	except: pass

	value = line[13:18].strip()
	try: json_file["Name"].append(str(value))
	except: pass

	value = line[18:24].strip()
	try: json_file["Chain"].append(str(value))
	except: pass

	value = line[24:36].strip()
	try: json_file["InterDG"].append(float(value))
	except: pass

	value = line[36:48].strip()
	try: json_file["InterDDG"].append(float(value))
	except: pass

	value = line[48:60].strip()
	try: json_file["NormTerDDG"].append(float(value))
	except: pass

	value = line[60:72].strip()
	try: json_file["IntraDG"].append(float(value))
	except: pass

	value = line[72:84].strip()
	try: json_file["IntraDDG"].append(float(value))
	except: pass

	value = line[84:96].strip()
	try: json_file["NormTraDDG"].append(float(value))
	except: pass

	value = line[96:107].strip()
	try: json_file["ChainAtoms"].append(int(value))
	except: pass

for element in json_file:
	print(len(json_file[element]))
json_file = json.dumps(json_file)

# Save as .json
file = sys.argv[1]
file = file.replace(".bals", ".json")
file = open(file, "w")
file.write(str(json_file))

# Save as .csv
# file = sys.argv[1]
# json_file = file.replace(".bals", ".json")
# print(json_file)
# df = pd.read_json(json_file)

# csv = df.to_csv()
# file = file.replace(".bals", ".csv")
# file = open(file, "w")
# file.write(csv.to_string())
