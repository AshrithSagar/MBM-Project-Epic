# encoding: utf-8

'''
Created on 21 Jul 2017

@author: Amaurys Ávila Ibarra
'''
import json, sys

from matplotlib.pyplot import show

from budeAlaScan.config import deltaG_key

from .plot_results import plot_mutations, plot_singles


def select_best_constellations(constellations, best_number):
    '''
    This will sort constellations by ddG and will return the best number of them.
    
    constellations -- dictionary containing the constellations.
    best_number -- Number of best constellations we want to select.
    
    return N best constellations sorted by ddG
    '''
    total_number = len(constellations)

    # If by any chance the number of best constellations required are equal or
    # greater than the total, return the constellations without doing anything.
    if total_number <= best_number:
        return constellations

    count_best = 0

    my_best_constellations = {}

    for key, value in sorted(constellations.items(), key=lambda kv: kv[1]):
        print(key, value[0])
        my_best_constellations[key] = value
#     for constellation in sorted(constellations.values()):
#         print(constellation[0])
#         print(constellation, constellations[constellation])
        count_best += 1
        if count_best >= best_number:
            break

    return my_best_constellations


def replot_ala_scan(fname):
    '''
    Re-plots the results from the alanine scanning app
    that were written into a json file.
    
    fname -- json file name 
    '''
    data_keys = ("pdb_id", "ala_scan", "mutants", "repack_scan")

    my_file = open(fname, 'r')

    try:
        my_data = json.load(my_file)
    except:
        print("Fatal Error: The file (%s) does not comply with the json format." % (fname))
        sys.exit(2)

    for key in data_keys:
        if key not in my_data:
            print("Fatal Error: The file (%s) does not seem to be produced by budeAlaScan." % (fname))
            sys.exit(2)

    pdb_id = my_data['pdb_id']
    ddGs = my_data['ala_scan']
    mutations_dGs = my_data['mutants']
    repack_ddGs = my_data['repack_scan']

#     mutations_dGs = select_best_constellations(mutations_dGs, 20)

    if len(ddGs) == 2:
        model_count = len(ddGs[1][deltaG_key])
    else:
        model_count = 1

    plot_singles(ddGs, pdb_id, model_count)
    wt_ddG = ddGs[0][deltaG_key]
    if mutations_dGs is not None:
        plot_mutations(repack_ddGs, mutations_dGs, pdb_id, model_count, wt_ddG)

    show()

    return
