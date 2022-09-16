#!/usr/bin/env python3
# encoding: utf-8

'''
Created on 5 Jul 2017

@author: Amaurys Ávila Ibarra
'''

from collections import OrderedDict
from copy import deepcopy
import itertools
from time import strftime, localtime

import isambard, ampal

from budeAlaScan import config
from budeAlaScan.config import deltaG_key, r_rNum_idx, r_ch_idx, r_rddG_idx

from . f_system import create_nested_dir


def group_residues_by_cutOff(selected_residues, my_ampal, cut_off, is_multimodel):
    """
    Get a list of selected residues and loop trough all them to group
    all the residues which are within the CA cut off distance from each other.
    
    If multiple models the CA distances will be consider from first model.
    
    selected_residues -- List of selected residues by ddG limit.
    my_ampal -- Ampal container or Assembly from PDB.
    cut_off -- Maximum distance among CA from residues within a constellation.
    
    Return a dictionary using the residues as key with 
    a list of all residues within the CA cut off limit. 
    """

    if is_multimodel:
        ref_ampal = my_ampal[0]
    else:
        ref_ampal = my_ampal

    my_groups = {}

    for res in selected_residues:
        r_ch_A = res[0]
        r_num_A = res[1:]
        my_groups[res] = [res]
        for res2 in selected_residues:
            if res == res2:
                continue
            r_ch_B = res2[0]
            r_num_B = res2[1:]
            ca_distance = isambard.ampal.geometry.distance(ref_ampal[r_ch_A][r_num_A]['CA'], ref_ampal[r_ch_B][r_num_B]['CA'])
            if  ca_distance <= cut_off:
                my_groups[res].append(res2)
    return my_groups


def create_auto_constellations(models_info, const_size, ddg_limit, my_ampal, cut_off, is_multimodel):
    '''
    Takes the information of single mutations for all residues and creates all
    possible constellations of const_size with residues which have ddG equal or greater than
    ddg_limit
    
    models_info -- list of residues' info
    const_size -- number of residues for each constellation
    ddg_limt -- Limit of ddG to contemplate a residue for constellation.
    my_ampal -- Ampal container or Assembly.
    cut_off -- Maximum distance among CA from residues within a constellation.
    '''

    selected_residues = []

    for residue in models_info:
        # Let's avoid the deltaG entry
        if residue is deltaG_key:
            continue
        if float(models_info[residue][r_rddG_idx]) >= ddg_limit:
            # Let's append the chain ID + residue number
            selected_residues.append(models_info[residue][r_ch_idx] + models_info[residue][r_rNum_idx])

    slcted_count = len(selected_residues)

    if slcted_count < 2:
        print("\nWarning: Not hot constellations could be made.")
        print("We could not find residues with ddG equal or greater to (%.2f) kJ/mol" % (ddg_limit))
        return None
    elif slcted_count < const_size:
        const_size = slcted_count
        print("\nWarning: We have reduced the constellation size to %d.\nNot enough residues were selected." % (slcted_count))
        print("If you want to select more residues.")
        print("Please decrease the BUDE ddG selection limit from (%.2f) kJ/mol." % (ddg_limit))

    my_groups = group_residues_by_cutOff(selected_residues, my_ampal, cut_off, is_multimodel)

    del selected_residues

    max_res_per_group = 0
    for a_group in my_groups:
        res_per_gruop = len(my_groups[a_group])
        if max_res_per_group < res_per_gruop:
            max_res_per_group = res_per_gruop

    if max_res_per_group == 1:
        print("\nWarning: No constellation was created.\nNot enough residues were selected.")
        print("If you want to select more residues.")
        print("Please increase the cut off distance selection limit from (%.2f) %s. " % (cut_off, config.angstrom_symb))
        print("and/or decrease the BUDE ddG selection limit from (%.2f) kJ/mol." % (ddg_limit))
        return None
    elif max_res_per_group < const_size:

        print("\nWarning: We have reduced the constellation size from '%s' to '%d'.\nNot enough residues were selected." % (const_size, max_res_per_group))
        print("If you want to select more residues.")
        print("Please increase the cut off distance selection limit from (%.2f) kJ/mol." % (cut_off))
        const_size = max_res_per_group

    constellations = []

    for a_group in my_groups:
        for a_const in itertools.combinations(sorted(my_groups[a_group]), const_size):
            constellations.append(a_const)

    constellations = set(constellations)

    return tuple(constellations)


def create_manual_constellations(mutants):
    '''
    This function will take the residue numbers given in the command line
    for each mutant and convert them into a list of constellations.
    The data is expected in the following format. R1,R2,Rn.
    Residues within a mutant is separated by comma and mutants
    are separated by space.
    
    mutants -- List of mutants 
    
    e.g:
    First mutan, residues 1, 2 and 3 chain A. Second, residues 7,8 chain B
    will create the constellations.
    constellations = [['A1','A2','A3'], ['B7',B8']] 
    
    Return a list of constellations
    '''
    constellations = []
    for mutan in mutants:
        a_constellation = tuple(mutan.split(','))
        constellations.append(a_constellation)

    return tuple(constellations)


def get_docking_unit_seq(unit_chains, my_ampal, is_multimodel):

    '''
    Get the a list of chains which form the docking unit, the ampal object and 
    return a dictionary with chain ids for keys and its values.
    The value is a list of two elements, the first is a flag set to false
    the second the sequence of the chain.
    
    unit_chains -- Chain(s) that form the docking unit.
    my_ampal -- ampal object from PDB.
    is_multimodel -- Whether we go a multiple model PDB or not. 
    
    Returns a dictionary sequences[Chain_ID] = [True, "SEQUENCE"]
    '''

    sequences = OrderedDict()

    for chain in unit_chains:
        # The first element is flag to change if that chain has a mutation or not.
        org_seq = [False]
        if is_multimodel:
            org_seq.append(my_ampal[0][chain].sequence)
        else:
            org_seq.append(my_ampal[chain].sequence)
        sequences[chain] = org_seq
    return sequences


def get_seqs_per_chain(my_ampal, is_multimodel, chain_ids=None):
    """
    Function for getting the sequence of an ampal object. If chain IDs are
    specified then it will return the sequeces for those chains, otherwise 
    it will return the sequence for the whole object.
    
    my_ampal -- Ampal container or Assembly.
    is_multimodel -- Multiple models or not.
    chain_id -- Chain IDs to get sequences from.
    
    returns sequence of selected chains
    """

    seq = OrderedDict()

    if is_multimodel:
        one_model = my_ampal[0]
    else:
        one_model = my_ampal

    if chain_ids:
        for an_id in chain_ids:
            seq[an_id] = one_model[an_id].sequence
    else:
        for a_chain  in one_model:
            seq[a_chain.id] = one_model[a_chain.id].sequence

    return seq


def get_dockingUnit_ampal(my_ampal, chain_ids, is_multimodel):
    """
    This function will take a ampal object from a pdb and return an ampal 
    object with the specified chains.
    
    my_ampal -- Ampal object from PDB, container or assembly.
    chain_ids -- Chains to select.
    
    return ampal object with selected chains.
    """

    if is_multimodel(my_ampal):
        docking_unit = ampal.AmpalContainer()
        for a_model in my_ampal:
            docking_unit_model = ampal.Assembly()
            for an_id in chain_ids:
                docking_unit_model.append(a_model[an_id])

            docking_unit.append(docking_unit_model)

    else:
        docking_unit = ampal.Assembly()
        for an_id in chain_ids:
            docking_unit.append(my_ampal[an_id])

    return docking_unit


def find_residue_position(res_number, sequence, monomer):
    '''
    Get a residue number, an ampal monomer and its sequence.
    It uses the length of the sequence to loop through the monomer
    until it finds the residue with residue number: res_number returning
    the position of the residue in the sequence.
    res_number:  Residue number to search for.
    sequence: String containing the sequence of the monomer, using one-letter residue naming. 
    monomer -- Ampal monomer 
    
    Return the position of the residue in the sequence if found,
    otherwise returns None.
    '''
    for pos in range(0, len(sequence)):
        if monomer[pos].id == res_number:
            return pos

    return None


def make_sequence_mutations(constellations, org_sequences, my_ampal, is_multimodel):
    '''
    Will create the mutated sequences for all hot constellations.
    It takes all constellations, parent sequences and ampal polymer object.
    It creates a dictionary for each constellation with its mutated sequence(s).
    
    constellations -- List of constellations to create mutated sequences for.
    org_sequences -- Original sequences for the polymer.
    my_ampal -- Ampal object to be mutated. It might be an Ampal Container or polymer.
    is_multimodel -- Whether we have a multiple mode pdb (container) or not (polymer).
    
    Returns a dictionary for each constellation. 
    mutations[mutation] = mutated_sequences
    '''
    mutations = {}

    flag_idx = 0
    seq_idx = 1

    if is_multimodel:
        my_polymer = my_ampal[0]
    else:
        my_polymer = my_ampal

    for constellation in constellations:

        mut_seq = deepcopy(org_sequences)

        for res in constellation:
            ch_id = res[0]
            res_number = res[1:]

            res_pos = find_residue_position(res_number, org_sequences[ch_id][seq_idx], my_polymer[ch_id])

            if res_pos is None:
                continue

            # We will modify sequence set flag to True.
            mut_seq[ch_id][flag_idx] = True

            # Change Residue to Alanine
            new_sequence = mut_seq[ch_id][seq_idx][:res_pos] + 'A' + mut_seq[ch_id][seq_idx][res_pos + 1 :]
            mut_seq[ch_id][seq_idx] = new_sequence

        mutations ["_".join(constellation)] = mut_seq

    return mutations


def remove_extra_constellations(constellations, pdb_id, ddg_limit, const_size):
    '''
    Print constellations that exceeds the Max number of constellation into a file.
    
    constellations -- Total of constellations created automatically.
    pdb_id -- Basename of the PDB file, PDB ID.
    ddg_limit -- delta delta G limit used for residue selection.
    const_size -- Size of the constellations, number of residues.
    
    Return only the first MaxAutoNumber constellations.  
    '''
    const_number = len(constellations)
    max_number = config.cfg.getint(config.general_sect, config.gen_max_auto_num_opt)
    time_stamp = strftime("%Y%m%d%H%M%S", localtime())

    extra_dir = config.cfg.get(config.dir_sect, config.dir_wrk_opt) + "/" + config.extra_const_dir

    fname = pdb_id + "_extraConstellations_" + time_stamp + ".txt"

    mssge = """
A total %d constellations were created with a size of %d residues.
The ddG limit used was %.4f.

We will do the first %d. 

The rest was written to the file 
[%s/%s]
in chunks of %d per line. If you still want to do that many
you could use the manual mode or re-run the auto mode after
increasing the (%s) in the [%s] section 
from the INI file [%s]

Please bear in mind if you decide for the second option
It will take a considerably longer time and the machine 
might crash.

Are you feeling lucky? :)

""" % (const_number, const_size, ddg_limit, max_number, config.extra_const_dir, fname, max_number,
       config.gen_max_auto_num_opt, config.general_sect, config.config_file)

    create_nested_dir(extra_dir)

    constellations_todo = constellations[:max_number]
    extra_constellations = constellations[max_number:]

    const_counter = 0

    full_fname = extra_dir + "/" + fname
    out_file = open(full_fname, 'w')

    out_file.write(mssge)
    print(mssge)

    for constellation in extra_constellations:
        out_file.write("%s " % (",".join(constellation)))
        const_counter += 1

        if const_counter == max_number:
            out_file.write("\n\n")
            const_counter = 0

    out_file.close()

    del extra_constellations

    return constellations_todo


def remove_extra_constellations_res(constellations, pdb_id, const_size):
    '''
    Print constellations that exceeds the Max number of constellation into a file.
    
    constellations -- Total of constellations created automatically.
    pdb_id -- Basename of the PDB file, PDB ID.
    const_size -- Size of the constellations, number of residues.
    
    Return only the first MaxAutoNumber constellations.  
    '''
    const_number = len(constellations)
    max_number = config.cfg.getint(config.general_sect, config.gen_max_auto_num_opt)
    time_stamp = strftime("%Y%m%d%H%M%S", localtime())

    extra_dir = config.cfg.get(config.dir_sect, config.dir_wrk_opt) + "/" + config.extra_const_dir

    fname = pdb_id + "_extraConstellations_" + time_stamp + ".txt"

    mssge = """
A total %d constellations were created with a size of %d residues.

We will do the first %d. 

The rest was written to the file 
[%s/%s]
in chunks of %d per line. If you still want to do that many
you could use the manual mode or re-run the auto mode after
increasing the (%s) in the [%s] section 
from the INI file [%s]

Please bear in mind if you decide for the second option
It will take a considerably longer time and the machine 
might crash.

Are you feeling lucky? :)

""" % (const_number, const_size, max_number, config.extra_const_dir, fname, max_number,
       config.gen_max_auto_num_opt, config.general_sect, config.config_file)

    create_nested_dir(extra_dir)

    constellations_todo = constellations[:max_number]
    extra_constellations = constellations[max_number:]

    const_counter = 0

    full_fname = extra_dir + "/" + fname
    out_file = open(full_fname, 'w')

    out_file.write(mssge)
    print(mssge)

    for constellation in extra_constellations:
        out_file.write("%s " % (",".join(constellation)))
        const_counter += 1

        if const_counter == max_number:
            out_file.write("\n\n")
            const_counter = 0

    out_file.close()

    return constellations_todo

