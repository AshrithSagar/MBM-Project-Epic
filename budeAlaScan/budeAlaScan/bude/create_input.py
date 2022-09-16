'''
Created on 14 Jun 2017

@author: Amaurys Ávila Ibarra

Module for all function to create the structure to run
BUDE Alanine scanning.
'''
from copy import deepcopy
import gc
import os

from isambard.modelling import scwrl

from budeAlaScan import config
from budeAlaScan.myutils.f_system import write_file, read_file
from budeAlaScan.myutils.misce import run_command, convert_aa


def form_docking_unit_name(chains_id, my_basename, mdl_padding_size, model_number=0, is_repacked=True):
    '''
    Creates the file name for a docking unit. 
    
    Name convention for docking units.
    Multi-model PDBID, Model Number Mn, Chain(s) ChAB and pdb extension.
    Single-model PDBID, Chain(s) ChAB and pdb extension.
    PDBID_Mn_Chs.pdb -> 1r8u_M01_ChA.pdb
    PDBID_Chs.pdb -> 1dan_ChHL.pdb
    
    chains_id -- Chain(s) that form the docking unit.
    my_basename -- Basename of the input file - PDB ID
    model_number -- Model number if multiple models PDB
    is_repacked -- Whether the was re-packed or not.
    
    Return the new name for the docking unit.
    '''

    if model_number > 0:
        if is_repacked:
            docking_unit_name = "{}_M{}_Ch{}_{}.pdb".format(my_basename, str(model_number).zfill(mdl_padding_size), chains_id, config.rpck_suffix)
        else:
            docking_unit_name = "{}_M{}_Ch{}.pdb".format(my_basename, str(model_number).zfill(mdl_padding_size), chains_id)
    else:
        if is_repacked:
            docking_unit_name = "{}_Ch{}_{}.pdb".format(my_basename, chains_id, config.rpck_suffix)
        else:
            docking_unit_name = "{}_Ch{}.pdb".format(my_basename, chains_id)

    return docking_unit_name


def form_mutated_unit_name(chains_id, my_basename, mutation, mdl_padding_size, model_number=0):
    '''
    Creates the file name for a mutated docking unit. 
    
    Name convention for docking units.
    Multi-model PDBID, Model Number Mn, Chain(s) ChAB, mutation and pdb extension.
    Single-model PDBID, Chain(s) ChAB, mutation and pdb extension.
    mutation -> A34_A36 
    PDBID_Mn_Chs_mutation.pdb -> 1r8u_M01_ChA_A34_A36.pdb
    PDBID_Chs_mutation.pdb -> 1dan_ChHL_H34_L36.pdb
    
    chains_id -- Chain(s) that form the docking unit.
    my_basename -- Basename of the input file - PDB ID
    mutation -- Residues that were mutated
    model_number -- Model number if multiple models PDB
    
    Return the new name for the mutated docking unit.
    '''

    if model_number > 0:
        docking_unit_name = "{}_M{}_Ch{}_{}.pdb".format(my_basename, str(model_number).zfill(mdl_padding_size), chains_id, mutation)
    else:
        docking_unit_name = "{}_Ch{}_{}.pdb".format(my_basename, chains_id, mutation)

    return docking_unit_name


def write_bude_pdb(pdb_name, pdb_content, ignore_hetatm=True):
    '''
    Write the PDB content to a file. The content of the PDB is a list of strings
    where each string is a line of the PDB.
    
    pdb_name -- File name for the PDB.
    pdb_content -- Content of the PDB.
    ignore_hetatm -- Whether we want to ignore HETATM or not.
    '''

    out_file = open(pdb_name, 'w')

    for line in pdb_content:
        if not line:
            continue

        if ignore_hetatm:
            if not line.startswith("HETATM"):
                out_file.write("%s\n" % line)
        else:
            out_file.write("%s\n" % line)

    out_file.write("END".ljust(80))
    out_file.close()


def get_docking_unit(model, chains):
    '''
    Get all chain(s) from a PDB model for a docking unit.
    
    model -- Ampal polymer object
    chains -- Chains that form the docking unit.
    
    Returns a list of lines that form all chains in the  docking unit.
    '''

    docking_unit = ""

    for ch in chains:
        docking_unit += model[ch].pdb

    return docking_unit.split('\n')


def print_docking_unit_pdb(my_model, unit_chains, pdb_basename, pdbs_dir, mdl_padding_size, model_number=0, ignore_hetatm=True, is_repacked=True):
    '''
    Prints a PDB for docking unit.
    
    my_model -- Ampal polymer object
    unit_chains -- Chains that for the docking unit.
    pdb_basename -- Basename of the input file - PDB ID
    pdbs_dir -- Directory where the file will be created.
    model_number -- Model number if multiple models PDB
    ignore_hetatm -- Whether we want ignore HETATM records or not.
    is_repacked -- Whether the polymer was re-packed or not.
    
    Returns the name of the PDB file for the docking unit.
    '''

    unit_chain_suffix = ''.join(unit_chains)
    unit_new_name = form_docking_unit_name(unit_chain_suffix, pdb_basename, mdl_padding_size, model_number, is_repacked=is_repacked)
    full_new_name = pdbs_dir + "/" + unit_new_name

    docking_unit = get_docking_unit(my_model, unit_chains)
    write_bude_pdb(full_new_name, docking_unit, ignore_hetatm=ignore_hetatm)

    return unit_new_name


def create_input_pdbs(my_ampal, is_multimodel, unit_chains, pdb_basename, pdbs_dir, is_repacked=True):
    '''
    Create PDBs for a BUDE docking unit with single or multiple models and
    create a list of file names for this docking unit.
    
    my_ampal -- Ampal Container or Polymer object
    is_multimodel -- Whether it is a multiple model unit or not.
    unit_chains -- Chains that form the docking unit.
    pdb_basename -- Basename of the input file - PDB ID
    pdbs_dir -- Directory where the file will be created.
    is_repacked -- Whether the Ampal object was re-packed or not.
    
    Returns a list with the created file names.
    '''

    ignore_hetatm = True
    doking_unit_new_names = []

    if is_multimodel:
        model_padding_size = len (str(len(my_ampal))) + 1
        model_number = 1
        for model in my_ampal:
            doking_unit_new_names.append(print_docking_unit_pdb(model, unit_chains, pdb_basename, pdbs_dir, model_padding_size, model_number, ignore_hetatm=ignore_hetatm, is_repacked=is_repacked))
            model_number += 1
    else:
        model = my_ampal
        model_padding_size = 0
        doking_unit_new_names.append(print_docking_unit_pdb(model, unit_chains, pdb_basename, pdbs_dir, model_padding_size, model_number=0, ignore_hetatm=ignore_hetatm, is_repacked=is_repacked))

    return doking_unit_new_names


def activate_seq_residues(seq_content, res_to_activate):
    """
    Fuction for activating all relevant resiudes in res_to_activate. 
    All other residues will deactivated.
    
    seq_content -- Sequence file content, list of lines.
    res_to_activate -- List of three-letter codes of residues to activate.
    
    returns a sequence file content with activated residues, list of lines.
    """

    updated_seq_file = []
    # Flags used by BUDE.
    active_flag = "true"
    inactive_flag = "false"

    for a_line in seq_content:
        if a_line.startswith(("%", "#", "END")):
            updated_seq_file.append(a_line)
            continue
        res_code = a_line.split()[0]

        if res_code in res_to_activate:
            update_line = "  {}            {}".format(res_code, active_flag)
        else:
            update_line = "  {}            {}".format(res_code, inactive_flag)
        updated_seq_file.append(update_line)

    return updated_seq_file


def  activate_seq_file(fname, res_to_activate):
    """
    This function will take a BUDE sequence file fname and activates 
    all residues in res_to_activate. It will write back the the activated 
    content to fname.
    
    fname -- BUDE sequence file name.
    res_to_activate -- List of three-letter codes of residues to activate.
     
    """

    seq_content = read_file(fname)
    activated_seq_file = activate_seq_residues(seq_content, res_to_activate)
    write_file(fname, activated_seq_file)
    return


def create_bude_sequence(full_pdb_fname, full_bseq_fname):
    """
    pdb_fname -- PDB file name to create BUDE sequence from.
    bseq_fname -- BUDE sequence file name.
    """
    seq_exe = config.cfg.get(config.exe_sect, config.exe_bseq_opt)
    bseq_cmd = "{} -i {} -o {}".format(seq_exe, full_pdb_fname, full_bseq_fname)
    run_command(seq_exe, bseq_cmd)
    if config.verbose:
        print(bseq_cmd)
    return


def create_input_bseq(pdb_fname, bseq_fname):
    """
    Create BUDE sequence for a pdb file.
    pdb_fname -- PDB file name to create BUDE sequence from.
    bseq_fname -- BUDE sequence file name.
    """

    pdb_full_fname = config.cfg.get(config.dir_sect, config.dir_l_src_opt) + "/" + pdb_fname

    create_bude_sequence(pdb_full_fname, bseq_fname)
    # create sequence command.
    seq_exe = config.cfg.get(config.exe_sect, config.exe_bseq_opt)
    bseq_cmd = "{} -i {} -o {}".format(seq_exe, pdb_full_fname, bseq_fname)
    run_command(seq_exe, bseq_cmd)
    if config.verbose:
        print(bseq_cmd)
    return


def create_dock_unit_list_file(files_names, pdb_basename, list_fname, files_dir):
    '''
    Create the file containing the list of file names for each docking unit.
    
    files_names -- List of file names.
    pdb_basename -- Base name of input pdb, (PDB ID).
    list_fname -- Name for list file.
    files_dir -- Directory where the files were created.
    '''
    os.chdir(config.cfg.get(config.dir_sect, config.dir_wrk_opt))
    scan_dir = config.cfg.get(config.dir_sect, config.dir_scan_opt)
    full_list_fname = scan_dir + "/" + pdb_basename + "_" + list_fname
    fnames = []
    for name in files_names:
        fnames.append(files_dir + "/" + name)
    write_file(full_list_fname, fnames)

    return


def create_dock_mut_units_list_fnames(mutations_srcs, unit_key):
    '''
    Create a single list adding mutation directory to file names for a docking unit.
  
    mutations_srcs -- List of file names per mutation.
    unit_key -- Docking unit we are interested in.
    '''

    unit_fnames = []

    for mutation in mutations_srcs[unit_key]:
        for fname in mutations_srcs[unit_key][mutation]:
            unit_fnames.append(mutation + "/" + fname)

    return unit_fnames


def get_mutated_sequence(mutated_sequences, my_model):
    '''
    This function will return the new mutated sequence for an 
    ampal polymer object. It takes the mutated sequence(s) and the polymer.
    
    The mutated sequences is a dictionary with the format below
    sequences[Chain_ID] = [True, "SEQUENCE"]
    
    mutated_sequences -- Sequences for the mutated chain(s).  
    my_model -- Ampal polymer object.
    
    Returns mutated sequence for the polymer.
    '''

    flag_idx = 0
    seq_idx = 1
    new_sequence = []

    for chain in my_model:
        if chain.id in mutated_sequences and mutated_sequences[chain.id][flag_idx]:
            new_sequence.append(mutated_sequences[chain.id][seq_idx])
        else:
            new_sequence.append(my_model[chain.id].sequence)
    return new_sequence


def get_mutations_sequences(mutations, my_ampal, is_multimodel):
    '''
    WE ARE NOT USING THIS FUNCTION AT THE MOMENT.
    
    Create all mutations sequences for all constellations.
    The mutations is a dictionary with the format below
    mutations[mutation] = sequences[Chain_ID] = [True, "SEQUENCE"]
    mutations -- Mutations with their monomers' mutated sequences.
    is_multimodel -- Whether we have a multiple model PDB or not.
    
    Returns a mutated sequences per constellation.
    '''

    mutations_sequences = {}

    if is_multimodel:
        my_model = my_ampal[0]
    else:
        my_model = my_ampal

    for mutation in mutations:
        mutations_sequences[mutation] = get_mutated_sequence(mutations[mutation], my_model)

    return mutations_sequences


def print_mutated_docking_unit_pdb(my_model, unit_chains, pdb_basename, pdbs_dir, mutation, model_padding_size, model_number=0, ignore_hetatm=True):
    '''
    Prints a PDB for docking unit.
    
    my_model -- Ampal polymer object
    unit_chains -- Chains that for the docking unit.
    pdb_basename -- Basename of the input file - PDB ID.
    pdbs_dir -- Directory where the file will be created.
    mutation -- mutation 
    model_number -- Model number if multiple models PDB.
    ignore_hetatm -- Whether we want ignore HETATM records or not.
    is_repacked -- Whether the polymer was re-packed or not.
    
    Returns the name of the PDB file for the docking unit.
    '''

    unit_chain_suffix = ''.join(unit_chains)

    unit_new_name = form_mutated_unit_name(unit_chain_suffix, pdb_basename, mutation, model_padding_size, model_number=model_number)
    full_new_name = pdbs_dir + "/" + unit_new_name

    docking_unit = get_docking_unit(my_model, unit_chains)
    write_bude_pdb(full_new_name, docking_unit, ignore_hetatm=ignore_hetatm)

    return unit_new_name


def create_mutations_pdbs(mutations, my_ampal, is_multimodel, rctr_chains, lig_chains, pdb_basename):
    '''
    Create PDBs for all mutated constellations
    The mutations is a dictionary with the format below
    mutations[mutation] = sequences[Chain_ID] = [True, "SEQUENCE"]
    my_ampal -- Ampal object that has been re-packed, it is either an AmpalContainer or polymer
    mutations -- Mutations with their monomers' mutated sequences.
    is_multimodel -- Whether we have a multiple model PDB or not.
    rctr_chains -- Chain(s) that form the receptor.
    lig_chains-- Chain(s) that form the ligand.
    pdb_basename -- Basename of the input file - PDB ID
    
    Return the two lists of names for the receptors and ligands.
    '''
    mut_srcs = {}
    receptor_new_names = {}
    ligand_new_names = {}

    pdbs_dir = config.cfg.get(config.dir_sect, config.dir_msrc_opt)
    ignore_hetatm = True

    os.chdir((config.cfg.get(config.dir_sect, config.dir_wrk_opt) + "/" + config.cfg.get(config.dir_sect, config.dir_scan_opt)))

    if is_multimodel:
        model_padding_size = len (str(len(my_ampal))) + 1
        org_model = deepcopy(my_ampal[0])
    else:
        model_padding_size = 0
        org_model = deepcopy(my_ampal)

    for mutation in mutations:

        receptor_new_names[mutation] = []
        ligand_new_names[mutation] = []

        mut_dir = pdbs_dir + "/" + mutation

        if is_multimodel:
            repacking_sequence = get_mutated_sequence(mutations[mutation], org_model)

            model_number = 1
            for model in my_ampal:

                repacked_model = scwrl.pack_side_chains_scwrl(model, repacking_sequence)

                receptor_new_names[mutation].append(print_mutated_docking_unit_pdb(repacked_model, rctr_chains, pdb_basename, mut_dir, mutation, model_padding_size, model_number, ignore_hetatm=ignore_hetatm))
                ligand_new_names[mutation].append(print_mutated_docking_unit_pdb(repacked_model, lig_chains, pdb_basename, mut_dir, mutation, model_padding_size, model_number, ignore_hetatm=ignore_hetatm))

                model_number += 1
                del repacked_model

        else:
            repacking_sequence = get_mutated_sequence(mutations[mutation], org_model)

            repacked_model = scwrl.pack_side_chains_scwrl(my_ampal, repacking_sequence)
            receptor_new_names[mutation].append(print_mutated_docking_unit_pdb(repacked_model, rctr_chains, pdb_basename, mut_dir, mutation, model_padding_size, model_number=0, ignore_hetatm=ignore_hetatm))
            ligand_new_names[mutation].append(print_mutated_docking_unit_pdb(repacked_model, lig_chains, pdb_basename, mut_dir, mutation, model_padding_size, model_number=0, ignore_hetatm=ignore_hetatm))
            del repacked_model

        gc.collect()

    mut_srcs.update({'receptors':receptor_new_names})
    mut_srcs.update({'ligands': ligand_new_names})

    del org_model

    return mut_srcs


def create_bude_seqs_constallations(constellations_srcs, mdl_number):
    """
    Create and activate BUDE sequence for all constellations' mutants.
    The active directory should the Scan directory. 'alaScan'
    
    constellations_srcs -- Dictionary with pdb mutant names per mutation.
    mdl_number -- number of models.
    """

    code3letter = [convert_aa(e) for e in config.res2activate]
    pdbs_dir = config.cfg.get(config.dir_sect, config.dir_msrc_opt)

    mdl_pad_size = len(str(mdl_number)) + 1
    mdl_suffix = "_M{}".format('1'.zfill(mdl_pad_size))

    for d_unit in constellations_srcs:
        dock_unit = constellations_srcs[d_unit]
        for mut in dock_unit:
            pdb_fname = dock_unit[mut][0]
            full_pdb_fname = "{}/{}/{}".format(pdbs_dir, mut, pdb_fname)

            seq_fname = pdb_fname.replace(mdl_suffix, "")
            seq_fname = seq_fname.replace(".pdb", ".bsql")
            full_seq_fname = "{}/{}/{}".format(pdbs_dir, mut, seq_fname)

            create_bude_sequence(full_pdb_fname, full_seq_fname)
            activate_seq_file(full_seq_fname, code3letter)

    return

####

