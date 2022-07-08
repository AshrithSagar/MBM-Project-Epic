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
	value = line[0:9].strip()
	try: json_file["Index"].append(int(value))
	except: pass

	value = line[10:14].strip()
	try: json_file["Number"].append(int(value))
	except: pass

	value = line[15:20].strip()
	try: json_file["Name"].append(str(value))
	except: pass

	value = line[21:25].strip()
	try: json_file["Chain"].append(str(value))
	except: pass

	value = line[26:39].strip()
	try: json_file["InterDG"].append(float(value))
	except: pass

	value = line[40:51].strip()
	try: json_file["InterDDG"].append(float(value))
	except: pass

	value = line[52:60].strip()
	try: json_file["NormTerDDG"].append(float(value))
	except: pass

	value = line[62:72].strip()
	try: json_file["IntraDG"].append(float(value))
	except: pass

	value = line[73:85].strip()
	try: json_file["IntraDDG"].append(float(value))
	except: pass

	value = line[86:101].strip()
	try: json_file["NormTraDDG"].append(float(value))
	except: pass

	value = line[102:107].strip()
	try: json_file["ChainAtoms"].append(int(value))
	except: pass

json_file = json.dumps(json_file)

# Save as .json
file = sys.argv[1]
file = file.replace(".bals", ".json")
file = open(file, "w")
file.write(str(json_file))

# Save as .csv
file = sys.argv[1]
json_file = file.replace(".bals", ".json")
print(json_file)
df = pd.read_json(json_file)
# csv = df.to_csv()
# file = file.replace(".bals", ".csv")
# file = open(file, "w")
# file.write(csv.to_string())
