#!/usr/bin/env python3
# encoding: utf-8

'''
Created on 12 Jun 2017

@author: Amaurys Ávila Ibarra
'''

import sys

from budeAlaScan import config
from budeAlaScan.myutils.misce import upper_my_list, create_unique_list, \
    expand_to_single_chars

from .. myutils.misce import get_legal_chains, is_multi_model

__all__ = ["check_input_common_errors_4modes", "check_input_errors_auto",
           "check_input_errors_manual", "check_input_errors_residues",
           "remove_duplicated_mutans"]


def check_docking_chains(legal_chains, docking_chains, pdb_name):
    '''
    Check that the chain ID(s) giving for a docking unit exist in the PDB.
    If the chain ID is not found it aborts execution.
    
    legal_chains -- List of chain(s) that are in the PDB
    docking_chains -- List of chain(s) for the docking unit.
    pdb_name -- PDB file name.
      
    '''
    for chain in docking_chains:
        if not is_chain_legal(legal_chains, chain, pdb_name):
            sys.stderr.write("FATAL ERROR: We could not find Chain ({}).\n".format(chain))
            sys.exit(2)

    return


def get_non_std_res_info(my_peptide):
    '''
    Get name and number from all non standard residues within a 
    peptide (chain)
    
    my_peptide -- Peptide object to check for non standard residues.
    
    returns a list of of residue(s) info.
    '''

    # residue index counter.
    res_idx = 0
    # list of none std residue indexes
    res_idxs = []
    # list of none std residue's info (name -> number)
    res_info = []

    for a_res in my_peptide.sequence:
        if a_res == 'X':
            res_idxs.append(res_idx)
        res_idx += 1

    for idx in res_idxs:
        an_info = "Name: {} -> No: {}".format(my_peptide[idx].mol_code, my_peptide[idx].id)
        res_info.append(an_info)

    return res_info


def non_std_in_model(model):
    '''
    Check if there is a non standard residue within a polymer.
    If any is found, returns True otherwise returns False.
    
    model -- Ampal polymer object.
    
    Returns true if none standard residue is found within the polymer, false otherwise.
    '''
    for chain in model:
        if 'X' in model[chain.id].sequence:
            sys.stderr.write("ERROR: We found non standard residue(s) in Chain ({}).\n".format(chain.id))
            sys.stderr.write("Please check the following residue(s):\n")
            for r_info in get_non_std_res_info(model[chain.id]):
                sys.stderr.write("{}\n".format(r_info))
            return True
    return False


def found_non_std_residue(my_ampal):
    '''
    Check if there is a non standard residue within an ampal container.
    If any is found, returns True otherwise returns False.
    
    my_ampal -- Ampal polymer or AmpalContainer object
    
    Returns true if none standard residue is found within the ampal object, false otherwise.
    '''

    if is_multi_model(my_ampal):
        for model in my_ampal:
            if non_std_in_model(model):
                return True

        return False
    else:
        if non_std_in_model(my_ampal):
            return True
        else:
            return False


def check_input_common_errors_4modes(pdb_name, receptor_chains, ligand_chains,
                                     my_ampal, res2activate):
    '''
    Check input data which is common for both modes, auto and manual.
    Checking that the ampal object from PDB is usable for ALA Scanning, at least two chains.
    Checking that the receptor and ligand chains are in the PDB.
    
    pdb_name -- PDB file name.
    receptor_chains -- Chain(s) that will form the receptor.
    ligand_chains -- Chain(s) that will form the ligand.
    my_ampal -- Ampal polymer or AmpalContainer object.
    res2activate -- One-letter code of residues to activate for rotamer correction.
    '''

    # if we want rotamer correction and residues to activate have been given
    # on the command line.
    if config.do_rotamer_correction and isinstance(res2activate, list):
        res2activate = create_unique_list(upper_my_list(expand_to_single_chars(res2activate)))
        illegal_aa_found = False
        for aa in res2activate:
            if aa not in config.legal_aa:
                illegal_aa_found = True
                sys.stderr.write("FATAL ERROR: Illegal one letter code found ({}).\n".format(aa))
        if illegal_aa_found:
            sys.stderr.write("Legal one letter code must be one of:\n[{}]\n".format(", ".join(config.legal_aa)))
            sys.exit(2)
        # if no error found on command line input, let's update the residues to
        # activate.
        config.res2activate = tuple(sorted(res2activate))

    # Let's check we have a sensible ampal object.
    if not is_ampal_object_ok(my_ampal):
        sys.stderr.write("FATAL ERROR: We could not make sense of the PDB file ({}).\n".format(pdb_name))
        sys.stderr.write("Is the file in PDB format? If so, does it have only one chain?\n")
        sys.exit(2)

    # Checking for non standard amino acids in the PDB.

    if found_non_std_residue(my_ampal):
        sys.stderr.write("FATAL ERROR: We found non standard residue(s) in PDB file ({}).\n".format(pdb_name))
        sys.exit(2)

    my_legal_chains = get_legal_chains(my_ampal)

    check_docking_chains(my_legal_chains, receptor_chains, pdb_name)
    check_docking_chains(my_legal_chains, ligand_chains, pdb_name)

    return


def check_input_errors_auto(deltaDelta_g, constellation_size):
    '''
    Check input data for the auto mode.
    Checking that the size of the constellations is at least of 2 residues.
    Checking that ddG is positive.
    
    deltaDelta_g -- ddG cut off limit.
    constellation_size -- Size of constellations, number of residues in each constellation.
    '''
    if deltaDelta_g <= 0:
        sys.stderr.write("FATAL ERROR:  ddG value must be greater than zero. ({:f}).\n".format(deltaDelta_g))
        sys.exit(2)

    if constellation_size < 2:
        sys.stderr.write("FATAL ERROR: Constellation size must be at least two residues. ({})\n".format(constellation_size))
        sys.exit(2)

    return


def remove_duplicated_mutans(mutants):
    '''
    Checking if our mutants have duplicates.
    If duplicates are found it will remove them.
    The equality follows the logic in the examples below
    since both cases will produce the same mutant.
    1st: A3,a4,A5 == A3,A4,A5
    2nd: A3,A4,A5 == A5,A3,A4
    
    mutants -- List of mutants [A3,A4,A5 ...].
     
    Returns a list of unique mutants A3,A4,A5 ...].
    '''

    unique_mutants = []

    for constellation in mutants:
        a_constellation = ",".join(sorted(constellation.upper().split(',')))
        if a_constellation not in unique_mutants:
            unique_mutants.append(a_constellation)

    # Let's give warning if duplicates are found.
    if len(unique_mutants) != len(mutants):
        print("Warning: We have deleted duplicated mutants.")
        print(mutants)
        print(unique_mutants)

    return unique_mutants


def remove_duplicated_residues(input_res):
    '''
    Checking if our residues have duplicates.
    If duplicates are found it will remove them.
    The equality follows the logic in the examples below
    case insensitive.
    1st: a4 == A4 
    
    input_res -- List of residues [A3, A4, A5, a5].
     
    Returns a list of unique residues in upper case [ A3,A4,A5].
    '''

    unique_residues = sorted(list(set(upper_my_list(input_res))))

    # Let's give warning if duplicates are found.
    if len(unique_residues) != len(input_res):
        print("Warning: We have deleted duplicated residues.")
        print(input_res)
        print(unique_residues)

    return unique_residues


def check_input_errors_residues(my_ampal, my_residues, ligand_chains, is_multimodel, constellation_size):

    '''
    Check input data for the manual mode.
    Checking that all chains from constellations are ligand chains.
    Checking all residue numbers from all constellations are present in 
    the ligand chains.
    A constellation is a list of residues ['A4', 'A5', 'A7']
    
    If a chain ID from one of the residues is not part of the declared ligand chain(s)
    of the residue number is not in the chain specified by chain ID, it will abort the
    execution.
    
    my_ampal -- Ampal polymer or AmpalContainer object
    my_residues -- List of residues 
    ligand_chains -- Chains that will form the ligand.
    is_multimodel -- Whether we have multiple model PDB or not.
    '''

    residues_number = len(my_residues)

    if constellation_size < 2 :
        sys.stderr.write("FATAL ERROR: Constellation size cannot be less than two.\n")
        sys.stderr.write("You might consider the 'scan' mode for single alanine mutants.\n")
        sys.exit(2)

    if residues_number < 2 :
        sys.stderr.write("FATAL ERROR: The number of residues cannot be less than two.\n")
        sys.stderr.write("You might consider the 'scan' mode for single alanine mutants.\n")
        sys.exit(2)

    error_found = False
    cap_chain = upper_my_list(ligand_chains)

    for residue in my_residues:

        chain_id = residue[0]

        if chain_id not in cap_chain:
            sys.stderr.write("Error: Chain [{}] from residue [{}] is not part of ligand chains(s) [{}].\n".format(chain_id, residue, ' '.join(ligand_chains)))
            error_found = True;
            break

        if not is_legal_residue_number(my_ampal, chain_id, residue[1:], is_multimodel):
            sys.stderr.write("Error: We could not found residue number ({}) in chain ({}).\n".format(residue[1:], chain_id))
            error_found = True
            break;

    if error_found:
        sys.stderr.write("FATAL ERROR: Residue not found in PDB or not part of ligand.\n")
        sys.exit(2)

    return


def check_input_errors_manual(my_ampal, constellations, ligand_chains, is_multimodel):
    '''
    Check input data for the manual mode.
    Checking that all chains from constellations are ligand chains.
    Checking all residue numbers from all constellations are present in 
    the ligand chains.
    A constellation is a list of residues ['A4', 'A5', 'A7']
    
    If a chain ID from one of the residues is not part of the declared ligand chain(s)
    of the residue number is not in the chain specified by chain ID, it will abort the
    execution.
    
    my_ampal -- Ampal polymer or AmpalContainer object
    constellations -- List constellations 
    ligand_chains -- Chains that will form the ligand.
    is_multimodel -- Whether we have multiple model PDB or not.
    '''

    cap_chain = upper_my_list(ligand_chains)

    for constellation in constellations:

        break_outer_loop = False

        for residue in constellation:

            chain_id = residue[0]

            if chain_id.upper() not in cap_chain:
                sys.stderr.write("Error: Chain [{}] from constellation [{}] is not part of ligand chains(s) [{}].\n".format(chain_id, ','.join(constellation), ' '.join(ligand_chains)))
                break_outer_loop = True;
                break

            if not is_legal_residue_number(my_ampal, chain_id.upper(), residue[1:], is_multimodel):
                sys.stderr.write("Error: We could not found residue ({}) in chain ({}).\n".format(residue[1:], chain_id))
                break_outer_loop = True
                break;

        if break_outer_loop:
            break

    if break_outer_loop:
        sys.stderr.write("FATAL ERROR: Residue from constellation not found in PDB or not part of ligand.\n")
        sys.exit(2)

    return


def is_ampal_object_ok(my_ampal):
    '''
    Checking that the ampal object is sensible for ALA scanning.
    For multi-model PDBs we check it is an ampal container and for both single and multi model
    PDBs we check that there are more than one molecule (chain within the PDB)
    
    my_ampal -- Ampal polymer or AmpalContainer object.
    
    Returns true if the ampal object is a polymer or AmpalContainer of polymers.
    '''

    ampal_ok = False

    # If multi-model check that we have an AmpalContainer
    if is_multi_model(my_ampal) and len(my_ampal[0]._molecules) > 1:
        ampal_ok = True

    else:
        # If single model check we have molecules.
        if len(my_ampal._molecules) > 1:
            ampal_ok = True

    return ampal_ok


def is_chain_legal(legal_chains, chain_id, pdb_filename):
    '''
    Checking that a chain ID is a legal chain
    
    legal_chains -- Chain IDs found in the PDB.
    chain_id -- Chain ID to check for.
    pdb_filename -- PDB file name
    
    Return true if chain_id is found in legal_chains, false otherwise.
    '''

    found_chain_id = True

    if chain_id.upper() not in legal_chains:
        sys.stderr.write("Error: Chain ID ({}) not found in PDB ({}).\n".format(chain_id, pdb_filename))
        sys.stderr.write("We found these chains in the PDB: [{}].\n".format(', '.join(legal_chains)))
        found_chain_id = False

    return found_chain_id


def is_legal_residue_number(my_ampal, chain_id, residue_number, is_multimodel):
    '''
    Check that residue is in the PDB with the specified chain.
    
    my_ampal -- Ampal polymer or AmpalContainer object.
    chain_id -- Chain ID the residue belong to.
    residue_number -- Residue number to find in the chain.
    is_multimodel -- Whether we have multiple model PDB or not.
    
    Return true if residue number is found within the chain id, false otherwise.
    '''

    found_it = False

    if is_multimodel:
        residues = my_ampal[0][chain_id]
    else:
        residues = my_ampal[chain_id]

    for residue in residues:
        if residue.id == residue_number:
            found_it = True
            break

    return found_it

