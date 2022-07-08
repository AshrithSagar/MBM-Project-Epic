"""
.bals to .csv converter
"""
import sys
import json
import csv

with open(sys.argv[1], 'r') as file:
    contents = file.readlines()

json = {
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
	try: json['Index'].append(int(value))
	except: pass

	value = line[10:14].strip()
	try: json['Number'].append(int(value))
	except: pass

	value = line[15:20].strip()
	try: json['Name'].append(str(value))
	except: pass

	value = line[21:25].strip()
	try: json['Chain'].append(str(value))
	except: pass

	value = line[26:39].strip()
	try: json['InterDG'].append(float(value))
	except: pass

	value = line[40:51].strip()
	try: json['InterDDG'].append(float(value))
	except: pass

	value = line[52:60].strip()
	try: json['NormTerDDG'].append(float(value))
	except: pass

	value = line[62:72].strip()
	try: json['IntraDG'].append(float(value))
	except: pass

	value = line[73:85].strip()
	try: json['IntraDDG'].append(float(value))
	except: pass

	value = line[86:101].strip()
	try: json['NormTraDDG'].append(float(value))
	except: pass

	value = line[102:107].strip()
	try: json['ChainAtoms'].append(int(value))
	except: pass

print(json)

file = sys.argv[1]
file = file.replace('.bals', '.json')
file = open(file, 'w')
file.write(str(json))

file = sys.argv[1]
file = file.replace('.bals', '.json')
file = open(file, 'w')
file.write(str(json))
