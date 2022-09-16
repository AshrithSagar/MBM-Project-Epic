#!/usr/bin/env python3
# encoding: utf-8

'''
Created on 14 Jul 2017

@author: aa6168
'''

import json
from os import chdir
from time import strftime, localtime

from numpy import subtract

from budeAlaScan import config
from budeAlaScan.config import deltaG_key, r_rNum_idx, r_rNam_idx, r_ch_idx, r_rddG_idx, r_sd_idx, r_sch_idx

from .bars import plot_model_bar, plot_mutations_bar
from .box import plot_models_box, plot_mutations_box


def __sort_by_key(key):
    """
    Simple function to sort ddGs dictionary by key. The keys of residues are integers.
    Turning them into a a padded string of padding_size leading zeros. Any other key is not
    relevant for our case.
    
    padding_size has to be larger than the number of digits of the residues number.
    10 should be enough that would mean a ten digits number.
    
    Return key if not integer, padded string if key is an integer.
    """
    padding_size = 10
    try:
        my_int_key = int(key)
    except ValueError:
        return key
    else:
        return '{0:0{1}}'.format(my_int_key, padding_size)


def do_bars_data(model):
    '''
    Process ddGs data from a model to create 
    tuples for the bars plot. The residue numbers are form
    with the chain ID + the number of the residue from the PDB
    ie. B33
    
    model is a dictionary of tuples where 'dG' key has the deltaG value
    and each residue has the an incremental index starting from '1'
    
    model -- Data from single mutants alanine scanning.
    
    Returns four tuples with residues' numbers, names, ddGs and standard deviations.
    '''
    r_numbers = ()
    r_names = ()
    r_ddGs = ()
    std_dev = ()
    sch_size = ()

    for key in sorted(model, key=__sort_by_key):
        if key == deltaG_key:
            continue
        r_numbers += (model[key][r_ch_idx] + model[key][r_rNum_idx],)
        r_names += (model[key][r_rNam_idx],)
        r_ddGs += (model[key][r_rddG_idx],)
        std_dev += (model[key][r_sd_idx],)
        sch_size += (model[key][r_sch_idx],)

    return r_numbers, r_names, r_ddGs, std_dev, sch_size


def do_box_data(model):
    '''
    Process ddGs data from a model to create 
    tuples for the box plot. The residue numbers are form
    with the chain ID + the number of the residue from the PDB
    ie. B33
    
    model is a dictionary of tuples where 'dG' key has the deltaG value
    and each residue has the an incremental index starting from '1'
    
    model -- Data from single mutants alanine scanning.
    
    Returns three tuples with residues' numbers, names, ddGs.
    '''
    r_numbers = ()
    r_names = ()
    r_ddGs = ()

    for key in sorted(model, key=__sort_by_key):
        if key == deltaG_key:
            continue
        r_numbers += (model[key][r_ch_idx] + model[key][r_rNum_idx],)
        r_names += (model[key][r_rNam_idx],)
        r_ddGs += (model[key][r_rddG_idx],)

    return r_numbers, r_names, r_ddGs


def plot_singles(ddGs, pdb_id, mdl_num=1):
    '''
    Plot single alanine mutants ddGs per residue.
    ddGs is a list with two elements if multiple models
    a single element if single model.
    The first element is the ddGs or average ddGs for each residue
    Average if multiple models
    ddGs = [{'dG':dG, 'N': residue}, {'dG':dG, 'N': models_residue}] 
    
    ddGs -- List with ddGs 
    pdb_id -- Basename of the input file, PDB ID
    mdl_num -- Number of models from the PDB
    '''

    is_avg = False

    if len(ddGs) > 1:
        is_avg = True

    r_numbers, r_names, r_ddGs, std_dev, sch_size = do_bars_data(ddGs[0])

    mol_dG = ddGs[0][deltaG_key]
    plot_model_bar(r_numbers, r_names, r_ddGs, std_dev, pdb_id, mol_dG, mdl_num, is_avg)

    if is_avg:
        r_numbers, r_names, r_ddGs = do_box_data(ddGs[1])
        plot_models_box(r_numbers, r_names, r_ddGs, pdb_id, mol_dG, mdl_num)

    return


def do_mutations_bars_data(mutations, mut_dGs, single_dG):

    # ddG = Mutant - WildType
    ddGs = ()

    for mutation in mutations:
        ddGs += (mut_dGs[mutation][0] - single_dG,)

    return ddGs


def do_mutations_box_data(mutations, mut_dGs, sin_models_dG):
    # ddG = Mutant - WildType
    ddGs = ()

    for mutation in mutations:
        ddGs += (tuple(subtract(mut_dGs[mutation][1], sin_models_dG)),)
    return ddGs


def plot_mutations(ddGs, mutations_dGs, pdb_id, mdl_cnt, wt_ddG):
    '''
    Plot hot constellations alanine mutants ddGs per mutation.
    ddGs is a list with two elements if multiple models
    a single element if single model.
    The first element is the ddGs or average ddGs for each residue
    Average if multiple models
    ddGs = [{'dG':dG, 'N': residue}, {'dG':dG, 'N': models_residue}] 
    
    The mutations dGs is a dictionary with that holds the dGs for each mutation.
    
    dGs = [dG, (dG_M1, dG_Mn)]
    mutations_dGs = {mutation: dGs} 
    
    ddGs -- List with ddGs 
    mutations_dGs -- deltaG values per mutation for single and multiple models PDBs.
    pdb_id -- Basename of the input file, PDB ID
    mdl_num -- Number of models from the PDB
    wt_ddG -- Wild Type ddG from the none repack scan.
    '''

    is_avg = False

    if len(ddGs) > 1:
        is_avg = True

    mutations = tuple(mutations_dGs)

    mutations_ddG = do_mutations_bars_data(mutations, mutations_dGs, wt_ddG)
    plot_mutations_bar(mutations, mutations_ddG, pdb_id, wt_ddG, mdl_cnt)

    if is_avg:
        mut_models_ddG = do_mutations_box_data(mutations, mutations_dGs, wt_ddG)
        plot_mutations_box(mutations, mut_models_ddG, pdb_id , wt_ddG, mdl_cnt)

    return


def write_ddGs(fout, ddGs, is_avg, models):
    """
    Write ddGs for the single mutants.
    
    fout -- file handler to write to.
    ddGs -- single mutants ddGs
    is_avg -- Whether we got an average ddGs or not (multiple models)
    models -- Header for the number of models.
    """

    r_numbers, r_names, r_ddGs, std_dev, sch_size = do_bars_data(ddGs[0])

    fout.write("%s,%s,%s,%s,%s\n" % ("ResNumber", "ResName", "ddGs", "std_dev", "sideCh_size"))
    for i in range(len(r_numbers)):
        fout.write("%s,%s,%.4f,%.4f,%i\n" % (r_numbers[i], r_names[i], r_ddGs[i], std_dev[i], sch_size[i]))

    if is_avg:
        r_numbers2, r_names2, r_ddGs2 = do_box_data(ddGs[1])
        fout.write("\n\n%s,%s,%s\n" % ("ResNumber", "ResName", ",".join(models)))
        for i in range(len(r_numbers2)):
            ddGs_str = []
            for ddG in r_ddGs2[i]:
                ddGs_str.append("%.4f" % ddG)
            fout.write("%s,%s,%s\n" % (r_numbers2[i], r_names2[i], ",".join(ddGs_str)))
    return


def write_plot_data(ddGs, mutations_dGs, repack_ddGs, pdb_id, unit_suffix,
                    model_count, args):
    '''
    Writes plot data and Alanine Scan Data to files.
    The plot data is written to a comma separated values file
    and the scan data is written into a json file.
    
    ddGs is a list with two elements if multiple models
    a single element if single model.
    The first element is the ddGs or average ddGs for each residue
    Average if multiple models
    ddGs = [{'dG':dG, 'N': residue}, {'dG':dG, 'N': models_residue}] 
    
    The mutations dGs is a dictionary with that holds the dGs for each mutation.
    
    dGs = [dG, (dG_M1, dG_Mn)]
    mutations_dGs = {mutation: dGs} 
    
    ddGs -- List with ddGs 
    mutations_dGs -- deltaG values per mutation for single and multiple models PDBs.
    pdb_id -- Basename of the input file, PDB ID
    unit_suffix -- Suffix for the docking unit, Rec or Lig
    model_count -- Number of models from the PDB
    args -- arguments that were passed in the command line.
    '''

    time_stamp = strftime("%Y%m%d%H%M%S", localtime())
    replot_full = config.cfg.get(config.dir_sect, config.dir_wrk_opt) + "/" + config.replot_dir
    chdir(replot_full)

    if args.scan_mode == 'auto':
        mode = "automatically"
    elif args.scan_mode == 'manual':
        mode = "manually"
    elif args.scan_mode == 'residues':
        mode = "given residues"
    elif args.scan_mode == 'scan':
        mode = "single mutants only"

    csv_fname = "{}{}_{}_{}_plot_data.csv".format(pdb_id, unit_suffix,
                                                  args.scan_mode, time_stamp)
    json_fname = "{}{}_{}_{}.json".format(pdb_id, unit_suffix, args.scan_mode,
                                          time_stamp)

    f_csv = open(csv_fname, 'w')

    f_json = open(json_fname, 'w')
    json_data = {"pdb_id": pdb_id, "ala_scan" : ddGs, "mutants": mutations_dGs,
                 "repack_scan": repack_ddGs}
    json.dump(json_data, f_json, indent=2)
    f_json.close()

    lig_chs = ", ".join(sorted ([x.upper() for x in args.ligand_chains]))
    rec_chs = ", ".join(sorted([x.upper() for x in args.receptor_chains]))

    comment = """# Hot constellations for PDB: %s
#
# PDB had %i model(s)
# Constellation selection was done %s.
# Receptor Chain(s) [%s]
# Ligand Chain[s) [%s]
#
#
""" % (pdb_id, model_count, mode, rec_chs, lig_chs)

    f_csv.write(comment)

    is_avg = False

    if len(ddGs) > 1:
        is_avg = True
        models = []
        for j in range (1, model_count + 1):
            models.append('M' + '{num:02d}'.format(num=j))
    else:
        models = None

    f_csv.write("Single mutants: \n")

    write_ddGs(f_csv, ddGs, is_avg, models)

    if repack_ddGs:
        f_csv.write("\n\nRepacked single mutants: \n")
        write_ddGs(f_csv, repack_ddGs, is_avg, models)

    if mutations_dGs is not None:
        mutations = tuple(mutations_dGs)
        mutations_ddG = do_mutations_bars_data(mutations, mutations_dGs, ddGs[0][deltaG_key])

        if is_avg:
            mut_models_ddG = do_mutations_box_data(mutations, mutations_dGs, ddGs[1][deltaG_key])

        f_csv.write("\n\n%s,%s\n" % ("Mutation", "ddG"))
        for k in range (len(mutations)):
            f_csv.write("%s,%.4f\n" % (mutations[k], mutations_ddG[k]))

        if is_avg:
            f_csv.write("\n\n%s,%s\n" % ("Mutation", ",".join(models)))
            for i in range (len(mutations)):
                ddGs_str = []
                for ddG in mut_models_ddG[i]:
                    ddGs_str.append("%.4f" % ddG)
                f_csv.write("%s,%s\n" % (mutations[i], ",".join(ddGs_str)))

    f_csv.close()

    if config.verbose:
        print("Plot data was written to: (%s)." % (csv_fname))
        print("Scan data was written to: (%s)." % (json_fname))

    return

