'''
Creating layout for Alanine Scanning in the local directory.

Created on 20 Jun 2017

@author: Amaurys Ávila Ibarra
'''

from os import path, chdir

from budeAlaScan import config
from budeAlaScan.bude import control
from budeAlaScan.bude.dock_libs import *

from . f_system import create_nested_dir, write_file


def create_layout():

    sources_full = config.cfg.get(config.dir_sect, config.dir_l_src_opt)
    libs_full = config.cfg.get(config.dir_sect, config.dir_l_blibs_opt)
    rslt_full = config.cfg.get(config.dir_sect, config.dir_l_rslt_opt)
    chim_full = config.cfg.get(config.dir_sect, config.dir_l_chi_opt)

    replot_full = config.cfg.get(config.dir_sect, config.dir_wrk_opt) + "/" + config.replot_dir

    if not path.isdir(sources_full):
        create_nested_dir(sources_full)
    if not path.isdir(libs_full):
        create_nested_dir(libs_full)
    if not path.isdir(rslt_full):
        create_nested_dir(rslt_full)
    if not path.isdir(chim_full):
        create_nested_dir(chim_full)
    if not path.isdir(replot_full):
        create_nested_dir(replot_full)

    if config.do_rotamer_correction:
        bude_seq_full = config.cfg.get(config.dir_sect, config.dir_l_seq_opt)
        if not path.isdir(bude_seq_full):
            create_nested_dir(bude_seq_full)

    ctrl_name = config.cfg.get(config.dir_sect, config.dir_scan_opt) + '/' + config.cfg.get(config.names_sect, config.names_b_ctrl_opt)
    zero_gen_name = libs_full + '/' + config.cfg.get(config.names_sect, config.names_b_gen_zero_opt)
    transf_name = libs_full + '/' + config.cfg.get(config.names_sect, config.names_b_transf_opt)
    ff_name = libs_full + '/' + config.cfg.get(config.names_sect, config.names_b_ff_opt)
    emc_name = libs_full + '/' + config.cfg.get(config.names_sect, config.names_b_bmc_opt)
    rotlib_name = libs_full + '/' + config.cfg.get(config.names_sect, config.names_b_rlib_opt)
    rotchoice_name = libs_full + '/' + config.cfg.get(config.names_sect, config.names_b_rchoice_opt)

    if not path.isfile(ctrl_name):
        write_file(ctrl_name, control.bude_ctrl)

    if not path.isfile(zero_gen_name):
        write_file(zero_gen_name, gen_zero.gen_zero)

    if not path.isfile(transf_name):
        write_file(transf_name, trans_rot.trans_rot)

    if not path.isfile(ff_name):
        write_file(ff_name, bhff.hybrid)

    if not path.isfile(emc_name):
        write_file(emc_name, bemc.bemc)

    if not path.isfile(rotlib_name):
        write_file(rotlib_name, rotamer_library.rotamer_library)

    if not path.isfile(rotchoice_name):
        write_file(rotchoice_name, rotamer_choice.rotamer_choice)
    return


def create_mutations_dirs(mutations):
    '''
    If there are no directories within the scan directory for mutations
    it will create two directories for mutations sources and results.
    
    It will also these options into the global configuration for the
    directories section if they do not exist.
    '''

    # Let's change dir to the scan directory
    chdir(config.cfg.get(config.dir_sect, config.dir_wrk_opt) + "/" + config.cfg.get(config.dir_sect, config.dir_scan_opt))

    muts_src = "mutations_src"
    muts_rslt = "mutations_rslt"

    for mutation in mutations:
        mut_src = muts_src + "/" + mutation
        mut_rslt = muts_rslt + "/" + mutation

        if not path.isdir(mut_src):
            create_nested_dir(mut_src)

        if not path.isdir(mut_rslt):
            create_nested_dir(mut_rslt)

    if not config.cfg.has_option(config.dir_sect, config.dir_msrc_opt):
        config.cfg.set(config.dir_sect, config.dir_msrc_opt, muts_src)

    if not config.cfg.has_option(config.dir_sect, config.dir_mrslt_opt):
        config.cfg.set(config.dir_sect, config.dir_mrslt_opt, muts_rslt)

    return
