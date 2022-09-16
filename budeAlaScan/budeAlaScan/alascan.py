#!/usr/bin/env python3
# encoding: utf-8

'''
Created on 10 Jul 2017

@author: Amaurys Ávila Ibarra
'''
import itertools
import os
import shutil
from time import strftime, localtime

from budeAlaScan import config
from budeAlaScan.bude.create_input import create_input_pdbs, \
create_dock_unit_list_file, create_mutations_pdbs, \
create_dock_mut_units_list_fnames, create_input_bseq, activate_seq_file, \
    create_bude_seqs_constallations
from budeAlaScan.bude.run_bude import run_bude, get_residues_ddg, \
get_constellations_delta_g
from budeAlaScan.errorparser.check_input_info import check_input_errors_auto, \
check_input_errors_manual, remove_duplicated_mutans, check_input_errors_residues, \
remove_duplicated_residues
from budeAlaScan.myutils.constellations import create_manual_constellations, \
create_auto_constellations, get_docking_unit_seq, remove_extra_constellations, \
make_sequence_mutations, remove_extra_constellations_res
from budeAlaScan.myutils.f_system import create_backup_fname
from budeAlaScan.myutils.layout import create_mutations_dirs
from budeAlaScan.myutils.misce import convert_aa
from budeAlaScan.myutils.repack import repackAmpal


def singles_scan(my_ampal, is_multimodel, my_basename, rctr_prefix, lig_prefix,
                 model_pad_size, receptor_chains, ligand_chains):
    """"
    It takes an ampal object, creates all single alanine mutants, calculates and parse ddG for all
    residues in the receptor and the ligand.
     
    my_ampal -- Ampal object it is either an AmpalContainer or Assembly
    is_multimodel -- True or False flag whether Ampal object is from multiple-model PDB or not.
    my_basename -- Basename of the input file - PDB ID
    rctr_prefix -- Prefix to add to the receptor's file names for SD and list.
    lig_prefix -- Prefix to add to the ligand's file names for SD and list.
    model_pad_size -- Size of padding for model number in file naming.
    receptor_chains -- Chain(s) that form the receptor docking unit. 
    ligand_chains -- -- Chain(s) that form the ligand docking unit. 
    
    
    Returns the ddGs for residues in the receptor and the ligand.
    """
    do_single_mutants(my_ampal, is_multimodel, my_basename, receptor_chains, ligand_chains)

    lig_ddgs = get_residues_ddg(my_basename, is_multimodel, ligand_chains, lig_prefix, model_pad_size)
    rec_ddgs = get_residues_ddg(my_basename, is_multimodel, receptor_chains, rctr_prefix, model_pad_size)

    return lig_ddgs, rec_ddgs


def residues_scan(my_ampal, is_multimodel, my_basename, rctr_prefix, lig_prefix,
                  args, model_pad_size, receptor_chains, ligand_chains):
    """"
    It takes an ampal object, creates all single alanine mutants, calculates and parse ddG for all
    residues in the receptor and the ligand. It will also create all possible constellations from 
    the given residues in the command line.
     
    my_ampal -- Ampal object it is either an AmpalContainer or Assembly
    is_multimodel -- True or False flag whether Ampal object is from multiple-model PDB or not.
    my_basename -- Basename of the input file - PDB ID
    rctr_prefix -- Prefix to add to the receptor's file names for SD and list.
    lig_prefix -- Prefix to add to the ligand's file names for SD and list.
    args -- Arguments given in the command line.
    model_pad_size -- Size of padding for model number in file naming.
    receptor_chains -- Chain(s) that form the receptor docking unit. 
    ligand_chains -- -- Chain(s) that form the ligand docking unit. 
    
    
    Returns the ddGs for residues in the receptor and the ligand as well as the ddG for the constellations.
    """

    constellations = []
    unique_residues = remove_duplicated_residues(args.my_residues)
    residues_number = len(unique_residues)

    check_input_errors_residues(my_ampal, unique_residues, args.ligand_chains, is_multimodel, args.constellation_size)

    if residues_number < args.constellation_size:
        print("WARNING: The number of residues given is less than the constellation size.")
        print("Constellation size will be change from {} to {}".format(args.constellation_size, residues_number))

        args.constellation_size = residues_number

    lig_ddgs, rec_ddgs = singles_scan(my_ampal, is_multimodel, my_basename, rctr_prefix, lig_prefix, model_pad_size, receptor_chains, ligand_chains)

    if config.there_is_scwrl:
        for a_const in itertools.combinations(unique_residues, args.constellation_size):
            constellations.append(a_const)

        max_number = config.cfg.getint(config.general_sect, config.gen_max_auto_num_opt)

        if len(constellations) > max_number:
            print("\nWOW: We over did it. Didn't we? :)")
            constellations = remove_extra_constellations_res(constellations, my_basename, args.constellation_size)
            print("Working on %d hot constellations ..." % (max_number))

        if constellations is not None:
            constellations_delta_g, repacked_lig_ddGs = do_constellations(my_ampal, is_multimodel, constellations, my_basename, receptor_chains, ligand_chains, lig_prefix, model_pad_size)
        else:
            constellations_delta_g = None
            repacked_lig_ddGs = None
    else:
        constellations_delta_g = None
        repacked_lig_ddGs = None
        print("\nWARNING: We are not doing constellations.")
        print("ISAMBARD could not find Scwrl.")

    del constellations, unique_residues
    return lig_ddgs, rec_ddgs, constellations_delta_g, repacked_lig_ddGs


def auto_scan(my_ampal, is_multimodel, my_basename, rctr_prefix, lig_prefix, args, model_pad_size, receptor_chains, ligand_chains):
    """"
    It takes an ampal object, creates all single alanine mutants, calculates and parse ddG for all
    residues in the receptor and the ligand. It will also create all possible constellations based
    on the ddG limit, constellation size and the residue's CA cut off distance.
     
    my_ampal -- Ampal object it is either an AmpalContainer or Assembly
    is_multimodel -- True or False flag whether Ampal object is from multiple-model PDB or not.
    my_basename -- Basename of the input file - PDB ID
    rctr_prefix -- Prefix to add to the receptor's file names for SD and list.
    lig_prefix -- Prefix to add to the ligand's file names for SD and list.
    args -- Arguments given in the command line.
    model_pad_size -- Size of padding for model number in file naming.
    receptor_chains -- Chain(s) that form the receptor docking unit. 
    ligand_chains -- -- Chain(s) that form the ligand docking unit. 
    
    
    Returns the ddGs for residues in the receptor and the ligand as well as the ddG for the constellations.
    """

    check_input_errors_auto(args.deltaDelta_g, args.constellation_size)
    lig_ddgs, rec_ddgs = singles_scan(my_ampal, is_multimodel, my_basename, rctr_prefix, lig_prefix, model_pad_size, receptor_chains, ligand_chains)

    if config.there_is_scwrl:
        constellations = create_auto_constellations(lig_ddgs[0], args.constellation_size, args.deltaDelta_g, my_ampal, args.cut_off, is_multimodel)

        max_number = config.cfg.getint(config.general_sect, config.gen_max_auto_num_opt)
        if (constellations is not None and len(constellations) > max_number):
            print("\nWOW: We over did it. Didn't we? :)")
            constellations = remove_extra_constellations(constellations, my_basename, args.deltaDelta_g, args.constellation_size)
            print("Working on %d hot constellations ..." % (max_number))

        if constellations is not None:
            constellations_delta_g, repacked_lig_ddGs = do_constellations(my_ampal, is_multimodel, constellations, my_basename, receptor_chains, ligand_chains, lig_prefix, model_pad_size)
        else:
            constellations_delta_g = None
            repacked_lig_ddGs = None

        del constellations
    else:
        constellations_delta_g = None
        repacked_lig_ddGs = None
        print("\nWARNING: We are not doing constellations.")
        print("ISAMBARD could not find Scwrl.")

    return lig_ddgs, rec_ddgs, constellations_delta_g, repacked_lig_ddGs


def manual_scan(my_ampal, is_multimodel, my_basename, rctr_prefix, lig_prefix, args, model_pad_size, receptor_chains, ligand_chains):
    """"
    It takes an ampal object, creates all single alanine mutants, calculates and parse ddG for all
    residues in the receptor and the ligand. It will also calculate the ddG for the given hot constellations.
     
    my_ampal -- Ampal object it is either an AmpalContainer or Assembly
    is_multimodel -- True or False flag whether Ampal object is from multiple-model PDB or not.
    my_basename -- Basename of the input file - PDB ID
    rctr_prefix -- Prefix to add to the receptor's file names for SD and list.
    lig_prefix -- Prefix to add to the ligand's file names for SD and list.
    args -- Arguments given in the command line.
    model_pad_size -- Size of padding for model number in file naming.
    receptor_chains -- Chain(s) that form the receptor docking unit. 
    ligand_chains -- -- Chain(s) that form the ligand docking unit. 
    
    
    Returns the ddGs for residues in the receptor and the ligand as well as the ddG for the constellations.
    """

    # Remove combinations from the command line which will produce the same mutant.
    mutants = remove_duplicated_mutans(args.my_constellations)
    constellations = create_manual_constellations(mutants)
    check_input_errors_manual(my_ampal, constellations, args.ligand_chains, is_multimodel)

    lig_ddgs, rec_ddgs = singles_scan(my_ampal, is_multimodel, my_basename, rctr_prefix, lig_prefix, model_pad_size, receptor_chains, ligand_chains)

    if config.there_is_scwrl:
        constellations_delta_g, repacked_lig_ddGs = do_constellations(my_ampal, is_multimodel, constellations, my_basename, receptor_chains, ligand_chains, lig_prefix, model_pad_size)
    else:
        constellations_delta_g = None
        repacked_lig_ddGs = None
        print("WARNING: We are not doing constellations.")
        print("isambard could not find Scwrl.")

    del constellations, mutants
    return lig_ddgs, rec_ddgs, constellations_delta_g, repacked_lig_ddGs


def do_single_mutants(my_ampal, is_multimodel, my_basename, receptor_chains,
                      ligand_chains, was_repack_done=False):
    """
    Creates all the input PDBs and runs BUDE calculating alanine scanning ddGs
    for single alanine mutants.
    
    my_ampal -- Ampal object it is either an AmpalContainer or Assembly
    is_multimodel -- True or False flag whether Ampal object is from multiple-model PDB or not.
    my_basename --  Basename of the input file - PDB ID
    receptor_chains -- Chain(s) that form the receptor docking unit. 
    ligand_chains -- Chain(s) that form the ligand docking unit.
    was_repack_done -- Flag to denote if there was repacking in the ampal object.
    """

    my_receptors = create_input_pdbs(my_ampal, is_multimodel, receptor_chains, my_basename, config.cfg.get(config.dir_sect, config.dir_l_src_opt), was_repack_done)
    my_ligands = create_input_pdbs(my_ampal, is_multimodel, ligand_chains, my_basename, config.cfg.get(config.dir_sect, config.dir_l_src_opt), was_repack_done)

    # We only need to create sequences the first time.
    if config.do_rotamer_correction and not was_repack_done:

        config.rec_bseq_fname = "{}_Ch{}.bsql".format(my_basename, "".join(receptor_chains))
        config.lig_bseq_fname = "{}_Ch{}.bsql".format(my_basename, "".join(ligand_chains))
        full_path = config.cfg.get(config.dir_sect, config.dir_l_seq_opt) + "/"
        rec_bseq_fname = full_path + config.rec_bseq_fname
        lig_bseq_fname = full_path + config.lig_bseq_fname

        create_input_bseq(my_receptors[0], rec_bseq_fname)
        create_input_bseq(my_ligands[0], lig_bseq_fname)

        code3letter = [convert_aa(e) for e in config.res2activate]
        activate_seq_file(rec_bseq_fname, code3letter)
        activate_seq_file(lig_bseq_fname, code3letter)

    create_dock_unit_list_file(my_receptors, my_basename, config.cfg.get(config.names_sect, config.names_rcpts_list_opt), config.cfg.get(config.dir_sect, config.dir_src_opt))
    create_dock_unit_list_file(my_ligands, my_basename, config.cfg.get(config.names_sect, config.names_ligs_list_opt), config.cfg.get(config.dir_sect, config.dir_src_opt))

    if len(my_receptors) > 1 and config.verbose:
        print("Starting single alanine mutants for %i models at %s." % (len(my_receptors), strftime(config.time_fmt, localtime())))

    # run_bude(receptors, ligands, src_dir=None, out_dir=None, rec_bseq_fname=None, lig_bseq_fname=None)
    if config.do_rotamer_correction:
        run_bude(receptors=my_receptors, ligands=my_ligands,
                 rec_bseq_fname=config.rec_bseq_fname,
                 lig_bseq_fname=config.lig_bseq_fname)
    else:
        run_bude(my_receptors, my_ligands)

    if len(my_receptors) > 1 and config.verbose:
        print("Finished single alanine mutants for %i models at %s." % (len(my_receptors), strftime(config.time_fmt, localtime())))
    elif config.verbose:
        print("We finished doing single alanine mutants at %s." % (strftime(config.time_fmt, localtime())))

    del my_receptors, my_ligands
    return


def do_constellations(my_ampal, is_multimodel, constellations, my_basename,
                      receptor_chains, ligand_chains, lig_prefix, model_pad_size):
    """
    Does the alanine scanning for multiple alanine mutants, constellations.
    
    my_ampal -- Ampal object it is either an AmpalContainer or Assembly
    is_multimodel -- True or False flag whether Ampal object is from multiple-model PDB or not.
    constellations -- constellations, multiple alanine mutants.
    my_basename -- Basename of the input file - PDB ID
    receptor_chains -- Chain(s) that form the receptor docking unit. 
    ligand_chains -- -- Chain(s) that form the ligand docking unit.
    lig_prefix -- Prefix to add to the ligand's file names for SD and list. 
    model_pad_size -- Size of padding for model number in file naming.
    
    Returns the ddGs the constellations.
    """

    # Let's change back to working directory.
    os.chdir(config.cfg.get(config.dir_sect, config.dir_wrk_opt))

    # Repacking ampal for ddG with mutants which will be repacked.
    repacked_ma = repackAmpal(my_ampal, is_multimodel)
    was_repack_done = True

    scan_dir = config.cfg.get(config.dir_sect, config.dir_scan_opt)

    # # Create backup for docking units file names list without repacking.
    rec_list_fname = my_basename + "_" + config.cfg.get(config.names_sect, config.names_rcpts_list_opt)
    full_list_fname = "{}/{}".format(scan_dir, rec_list_fname)
    bck_fname = "{}/{}".format(scan_dir, create_backup_fname(rec_list_fname, scan_dir))
    # Let's backup the set of docking unit names without repacking
    shutil.copyfile(full_list_fname, bck_fname)

    lig_list_fname = my_basename + "_" + config.cfg.get(config.names_sect, config.names_ligs_list_opt)
    full_list_fname = "{}/{}".format(scan_dir, lig_list_fname)
    bck_fname = "{}/{}".format(scan_dir, create_backup_fname(lig_list_fname, scan_dir))
    # Let's backup the set of docking unit names without repacking
    shutil.copyfile(full_list_fname, bck_fname)

    # Let's run the single mutants again for the repack ampal to do ddG with mutation.
    do_single_mutants(repacked_ma, is_multimodel, my_basename, receptor_chains,
                      ligand_chains, was_repack_done)

    repacked_lig_ddgs = get_residues_ddg(my_basename, is_multimodel, ligand_chains,
                                lig_prefix, model_pad_size, was_repack_done)

    lig_sequences = get_docking_unit_seq(ligand_chains, repacked_ma, is_multimodel)

    mutations_seqs = make_sequence_mutations(constellations, lig_sequences, repacked_ma, is_multimodel)

    my_mutations = tuple(mutations_seqs)

    if config.verbose:
        print("Starting %i Hot Constellations at %s." % (len(my_mutations), strftime(config.time_fmt, localtime())))

    create_mutations_dirs(my_mutations)
    constellations_srcs = create_mutations_pdbs(mutations_seqs, repacked_ma, is_multimodel, receptor_chains, ligand_chains, my_basename)

    constellations_receptors = create_dock_mut_units_list_fnames(constellations_srcs, 'receptors')
    constellations_ligands = create_dock_mut_units_list_fnames(constellations_srcs, 'ligands')

    # Let's create & activate BUDE sequences.
    if config.do_rotamer_correction:
        create_bude_seqs_constallations(constellations_srcs, len(my_ampal))

    create_dock_unit_list_file(constellations_receptors, my_basename + "_Mutations", config.cfg.get(config.names_sect, config.names_ligs_list_opt), config.cfg.get(config.dir_sect, config.dir_msrc_opt))
    create_dock_unit_list_file(constellations_ligands, my_basename + "_Mutations", config.cfg.get(config.names_sect, config.names_rcpts_list_opt), config.cfg.get(config.dir_sect, config.dir_msrc_opt))

    for mutation in my_mutations:

        if config.verbose:
            print("Doing Hot Constellation: %s" % (mutation))

        src_dir = config.cfg.get(config.dir_sect, config.dir_msrc_opt) + "/" + mutation
        dest_dir = config.cfg.get(config.dir_sect, config.dir_mrslt_opt) + "/" + mutation + "/"
        mdl_pad_size = len(str(len(my_ampal))) + 1
        mdl_suffix = "_M{}".format('1'.zfill(mdl_pad_size))

        rec_bseq_fname = constellations_srcs['receptors'][mutation][0].replace(mdl_suffix, "")
        rec_bseq_fname = rec_bseq_fname.replace(".pdb", ".bsql")
        # full_rec_bseq_fname = "{}/{}".format(src_dir, rec_bseq_fname)

        lig_bseq_fname = constellations_srcs['ligands'][mutation][0].replace(mdl_suffix, "")
        lig_bseq_fname = lig_bseq_fname.replace(".pdb", ".bsql")
        # full_lig_bseq_fname = "{}/{}".format(src_dir, rec_bseq_fname)

        if config.do_rotamer_correction:
            run_bude(constellations_srcs['receptors'][mutation],
                     constellations_srcs['ligands'][mutation],
                     src_dir, dest_dir, rec_bseq_fname, lig_bseq_fname, src_dir)
        else:
            run_bude(constellations_srcs['receptors'][mutation],
                     constellations_srcs['ligands'][mutation],
                     src_dir, dest_dir)

    constellations_delta_g = get_constellations_delta_g(my_mutations, my_basename, is_multimodel, ligand_chains, lig_prefix, model_pad_size)

    del constellations_ligands, constellations_receptors
    del my_mutations, lig_sequences, mutations_seqs, constellations_srcs

    if config.verbose:
        print("Finished Hot Constellations at %s." % (strftime(config.time_fmt, localtime())))

    return constellations_delta_g, repacked_lig_ddgs

