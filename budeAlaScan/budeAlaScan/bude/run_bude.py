'''
Created on 30 Jun 2017

@author: Amaurys Ávila Ibarra
'''

import os, sys, shutil, glob

from budeAlaScan import config
from budeAlaScan.config import b_rddG_idx, b_sch_idx, b_sd_idx, r_rddG_idx
from budeAlaScan.config import deltaG_key, b_idx, b_rNum_idx, b_rNam_idx, b_ch_idx
from budeAlaScan.myutils.f_system import write_file, read_file
from budeAlaScan.myutils.misce import run_command, print_std_error

# Tag to find the deltaG result in the BUDE output.
# used by: parse_mutation_results, parse_models_ddg, parse_getsd_results
__delta_g_tag = "# WT InterDG: "


def parse_mutation_results(results):
    '''
    Parse results from BUDE alaScan for a model
    of a mutant. We are only interested in the deltaG of the interaction.
    
    results -- is a list of lines from BUDE alaScan.
    
    returns the dG of the docking unit.
    '''

    for line in results:
        if line.startswith(__delta_g_tag):
            delta_g = float(line[len(__delta_g_tag):])
            break

    return delta_g


def parse_mutation_models(filenames):
    '''
    Parse mutation models to get the deltaG.
    
    filenames -- List of file names for each model BUDE alanine scan result.
    
    returns a tuple with each deltaG per model.
    '''

    models_deltaG = ()

    for filename in filenames:
        f_content = read_file(filename)

        for line in f_content:
            if line.startswith(__delta_g_tag):
                models_deltaG += (float(line[len(__delta_g_tag):]),)
                break

    return models_deltaG


def parse_models_ddg(filenames):
    '''
    Parse the ddG for multiple models it creates a dictionary 
    
    filenames -- List of file names for each model BUDE alanine scan result.
    returns a dictionary from 1 to N number of residues plus dG.
    
    residue = (residue_number, residue_name, chain_ID, (dG_M1, dG_Mn), sidechain_size) 
    returns models_info = {'dG':dG, 'N': residue}
    
    '''

    my_models = {}

    for filename in filenames:
        f_content = read_file(filename)

        for line in f_content:
            if line.startswith(__delta_g_tag):
                delta_g = float(line[len(__delta_g_tag):])
                if deltaG_key in my_models:
                    my_models [deltaG_key] = my_models[deltaG_key] + (float(delta_g),)
                else:
                    my_models [deltaG_key] = (float(delta_g),)
                continue

            if line.startswith('#'):
                continue

            res = line.split()

            if res[b_idx] in my_models:
                my_models[res[b_idx]] = (res[b_rNum_idx], res[b_rNam_idx],
                                         res[b_ch_idx],
                                         my_models[res[b_idx]][r_rddG_idx] +
                                         (float(res[b_rddG_idx]),),
                                         int(res[b_sch_idx]))
            else:
                my_models[res[b_idx]] = (res[b_rNum_idx], res[b_rNam_idx],
                                         res[b_ch_idx],
                                         (float(res[b_rddG_idx]),),
                                         int(res[b_sch_idx]))

    return my_models


def parse_getsd_results(results, is_multimodel):
    '''
    Parse results from the getSD utility or BUDE alaScan.
    results -- is a list of lines from the getSD utility or the BUDE alaScan
    is_multimodel -- whether we have multiple model PDB or not.
    
    It returns a dictionary with two elements, the first is the dG of the interaction and 
    the second is a list of lists holding info for all residues.
    
    
    Each nested list has the following indexes.
    0: Residue Number
    1: Residue Name
    2: Chain
    3: ddG
    4: sidechain size
    5: SD, it will be 0.0 if single model.
    
    returns molecule_info = {'delta_g' : delta_g, 'residues' : residues[ residue[0, 1, 2, 3, 4, 5] ] }
    '''
    single_mutation_results = {}

    for line in results:
        if line.startswith(__delta_g_tag):
            single_mutation_results [deltaG_key] = float(line[len(__delta_g_tag):])
            continue
        # Let's avoid empty or comment lines.
        if not line or line.startswith('#'):
            continue

        my_fields = line.split()
        # If no a multi model there is no SD.
        if not is_multimodel:
            my_fields.append(0.0)

        single_mutation_results[my_fields[b_idx]] = (my_fields[b_rNum_idx],
                                                     my_fields[b_rNam_idx],
                                                     my_fields[b_ch_idx],
                                                     float(my_fields[b_rddG_idx]
                                                           ),
                                                     int(my_fields[b_sch_idx]),
                                                     float(my_fields[b_sd_idx]))

    return single_mutation_results


def run_bude(receptors, ligands, src_dir=None, out_dir=None,
             rec_bseq_fname=None, lig_bseq_fname=None, bseq_dir=None):
    '''
    Run BUDE for each receptor-ligand pair to calculate ddGs.
    If the source and output directories are None
    We are doing single mutants and the directory for the source and output
    will be those from the single mutants. Otherwise we use src_dir and out_dir
    to passed to BUDE in the command line.
    
    receptors -- List of receptors' file names.
    ligands -- List of ligands' file names.
    src_dir -- Directory containing PDBs for receptors and ligands.
    out_dir -- Directory to output BUDE's results.
    rec_bseq_fname -- BUDE sequence file name for the receptor.
    lig_bseq_fname -- BUDE sequence file name for the ligand.
    bseq_dir -- directory where to find the BUDE sequences.
    '''

    my_exe = config.cfg.get(config.exe_sect, config.exe_bude_opt)
    ctrl_file = config.cfg.get(config.names_sect, config.names_b_ctrl_opt)

    if src_dir is None:
        src_dir = config.cfg.get(config.dir_sect, config.dir_src_opt)
    if out_dir is None:
        out_dir = config.cfg.get(config.dir_sect, config.dir_rslt_opt)
    if bseq_dir is None:
        bseq_dir = config.cfg.get(config.dir_sect, config.dir_seq_opt)

    chim_dir = config.cfg.get(config.dir_sect, config.dir_chi_opt)

    os.chdir(config.cfg.get(config.dir_sect, config.dir_wrk_opt) + "/" + config.cfg.get(config.dir_sect, config.dir_scan_opt))

    for receptor, ligand in zip(receptors, ligands):
        if config.do_rotamer_correction and rec_bseq_fname and lig_bseq_fname:
            bude_cmd = ("{0} -f {1} -P {2}/{3} -I {2}/{4} -R {5} "
            "-S {6}/{7} -Q {6}/{8}".format(my_exe, ctrl_file, src_dir, receptor, ligand,
                                out_dir, bseq_dir, rec_bseq_fname,
                                lig_bseq_fname)
            )
        else:
            bude_cmd = ("{0} -f {1} -P {2}/{3} -I {2}/{4} -R {5} "
            "-a 00".format(my_exe, ctrl_file, src_dir, receptor, ligand,
                                 out_dir)
            )
        run_command(my_exe, bude_cmd)
        if config.verbose:
            print(bude_cmd)
    # Let' move all chimera scripts.
    for file in glob.glob("*.com"):
        dest = "%s/%s" % (chim_dir, file)
        shutil.move(file, dest)

    return


def get_results_list(r_pattern):
    '''
    Creates a list of files name that matches the require pattern
    If not files are matched, it aborts execution.
    
    r_pattern -- pattern to match for results files.
    
    Return a list with files names.
    '''
    rslt_list = glob.glob(r_pattern)

    if len(rslt_list) == 0:
        my_error = "Fatal Error: We could not find any results with this pattern (%s).\n" % (r_pattern)
        my_error += "Searched path: [%s]." % (os.getcwd())
        print_std_error(my_error)
        sys.exit(4)
    return rslt_list


def get_residues_ddg(pdb_basename, is_multimodel, unit_chains, unit_prefix,
                     model_pad_size, was_repack_done=False):
    '''
    We will execute the getSD utility if we have a multiple model PDB
    otherwise we will use BUDE ALA scanning results.
    
    For multiple models PDBs it also returns the information for each model.
    
    The results are parsed and a list of lists with the results for each residue is returned.
    pdb_basename -- base name of the PDB file given to scan.
    is_multimodel -- whether we have multiple model PDB or not.
    unit_chains -- list of chain(s) that form the docking unit.
    unit_prefix -- Prefix to add to the file names of the list of results and standard deviation.
    was_repack_done -- Flag to denote if there was repacking in the ampal object.
    
    returns a two elements list of dictionaries from 1 to N number of residues plus dG
    the first element is the average or single model ddGs results
    the second element is the results per model.
    single_results = [{'dG':dG, 'N': residue}, {'dG':dG, 'N': models_residue}] 
    '''

    rslt_dir = config.cfg.get(config.dir_sect, config.dir_l_rslt_opt)
    work_dir = config.cfg.get(config.dir_sect, config.dir_wrk_opt)
    os.chdir((work_dir + "/" + rslt_dir))

    my_exe = config.cfg.get(config.exe_sect, config.exe_getsd_opt)

    single_mut_results = []

    if is_multimodel:
        match_padding = '?' * model_pad_size

        sd_fname = "%s_%s_Ch%s_SD.bals" % (unit_prefix, pdb_basename, ''.join(unit_chains))
        list_fname = "%s_%s_Ch%s_results.list" % (unit_prefix, pdb_basename, ''.join(unit_chains))

        if was_repack_done:
            rslt_pattern = "{}_M{}_Ch{}_{}*_{}0*.bals".format(pdb_basename,
                                                              match_padding,
                                                              ''.join(unit_chains),
                                                              config.rpck_suffix,
                                                              unit_prefix[0])
        else:
            rslt_pattern = "{}_M{}_Ch{}*_{}0*.bals".format(pdb_basename,
                                                           match_padding,
                                                           ''.join(unit_chains),
                                                           unit_prefix[0])

        rslt_list = get_results_list(rslt_pattern)
        write_file(list_fname, rslt_list)

        sd_cmd = "%s %s" % (my_exe, list_fname)
        sd_proc = run_command(my_exe, sd_cmd)
        write_file(sd_fname, sd_proc.stdout)

        single_mut_results.append(parse_getsd_results(sd_proc.stdout.split('\n'), is_multimodel))

        single_mut_results.append(parse_models_ddg(rslt_list))

    else:
        if was_repack_done:
            rslt_pattern = "{}_Ch{}_{}*_{}0*.bals".format(pdb_basename,
                                                      ''.join(unit_chains),
                                                      config.rpck_suffix,
                                                      unit_prefix[0])
        else:
            rslt_pattern = "{}_Ch{}_*_{}0*.bals".format(pdb_basename,
                                                    ''.join(unit_chains),
                                                    unit_prefix[0])
        rslt_list = get_results_list(rslt_pattern)
        single_mut_results.append(parse_getsd_results(read_file(rslt_list[0]), is_multimodel))

    return single_mut_results


def get_constellations_delta_g(mutations, pdb_basename, is_multimodel, unit_chains, unit_prefix, model_pad_size):
    '''
    For each mutation (constellation).
    We will execute the getSD utility if we have a multiple model PDB
    otherwise we will use BUDE ALA scanning results.
    
    The results are parsed and a dictionary is create with an entry per mutation
    Each entry will have list of two elements. The first element will be a float
    deltaG for single-model PDBs and average deltaG for multiple-models PDBs.
    The second element is present only for multiple models PDBs.
    It will be a tuple with the deltaG for each model. 
        
    mutations -- List of mutations.
    pdb_basename -- base name of the PDB file given to scan.
    is_multimodel -- whether we have multiple model PDB or not.
    unit_chains -- list of chain(s) that form the docking unit.
    unit_prefix -- Prefix to add to the file names of the list of results and standard deviation.
    
    returns a dictionary of mutations and their respective dGs
    '''

    rslt_dir = config.cfg.get(config.dir_sect, config.dir_scan_opt) + "/" + config.cfg.get(config.dir_sect, config.dir_mrslt_opt)
    work_dir = config.cfg.get(config.dir_sect, config.dir_wrk_opt)
    os.chdir((work_dir + "/" + rslt_dir))

    my_exe = config.cfg.get(config.exe_sect, config.exe_getsd_opt)
    constellations_delta_g = {}
    for mutation in mutations:
        models_delta_g = []
        if is_multimodel:
            match_padding = '?' * model_pad_size
            sd_fname = "%s_%s_Ch%s_%s_SD.bals" % (unit_prefix, pdb_basename, ''.join(unit_chains), mutation)
            list_fname = "%s_%s_Ch%s_%s_results.list" % (unit_prefix, pdb_basename, ''.join(unit_chains), mutation)
            rslt_pattern = "{}/{}_M{}_Ch{}_{}*.bals".format(mutation, pdb_basename, match_padding, ''.join(unit_chains), mutation)

            rslt_list = get_results_list(rslt_pattern)
            write_file(list_fname, rslt_list)

            sd_cmd = "%s %s" % (my_exe, list_fname)
            sd_proc = run_command(my_exe, sd_cmd)
            write_file(sd_fname, sd_proc.stdout)

            models_delta_g.append(parse_mutation_results(sd_proc.stdout.split('\n')))

            models_delta_g.append(parse_mutation_models(rslt_list))

        else:
            rslt_pattern = "%s/%s_Ch%s_%s*.bals" % (mutation, pdb_basename, ''.join(unit_chains), mutation)
            rslt_list = get_results_list(rslt_pattern)
            models_delta_g.append(parse_mutation_results(read_file(rslt_list[0])))

        constellations_delta_g[mutation] = models_delta_g

    return constellations_delta_g

