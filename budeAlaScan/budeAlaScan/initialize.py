'''
This module will read

Created on 19 Jun 2017

@author: Amaurys Ávila Ibarra
'''

from configparser import ConfigParser, ExtendedInterpolation
import os, errno, stat
from shutil import copy2

__config_file = "%s/.budeAlaScan.ini" % os.path.expanduser("~")

# Ini file sections
__exe_sect = 'Executables'
__dir_sect = 'Directories'
__names_sect = 'File Names'
__general_sect = 'General'

# Ini file General options.
__gen_max_auto_num_opt = "MaxAutoNumber"
__gen_max_auto_num_val = '20'

# ini file dir options
__dir_wrk_opt = "Work"
__dir_wrk_val = "%s/budeAlaScan" % (os.path.expanduser("~"))
__dir_scan_opt = "Scan"
__dir_scan_val = "alaScan"

__dir_blibs_opt = "BudeLibs"
__dir_blibs_val = "libs"
__dir_rslt_opt = "BudeResults"
__dir_rslt_val = "results"

__dir_src_opt = "PdbSources"
__dir_src_val = "sources"

__dir_seq_opt = "BudeSeqs"
__dir_seq_val = "budeSeqs"

__dir_chi_opt = "Chimera"
__dir_chi_val = "chimeraScripts"

__dir_l_src_opt = "LongPdbSources"
__dir_l_src_val = "${%s}/${%s}" % (__dir_scan_opt, __dir_src_opt)

__dir_l_seq_opt = "LongBudeSeqs"
__dir_l_seq_val = "${%s}/${%s}" % (__dir_scan_opt, __dir_seq_opt)

__dir_l_blibs_opt = "LongBudeLibs"
__dir_l_blibs_val = "${%s}/${%s}" % (__dir_scan_opt, __dir_blibs_opt)

__dir_l_rslt_opt = "LongResults"
__dir_l_rslt_val = "${%s}/${%s}" % (__dir_scan_opt, __dir_rslt_opt)

__dir_l_chi_opt = "LongChimera"
__dir_l_chi_val = "${%s}/${%s}" % (__dir_scan_opt, __dir_chi_opt)

__dir_cpp_opt = "CppCode"
__dir_cpp_val = "cppCode"

__dir_bude_opt = "BUDE"
__dir_bude_val = "${%s}/%s/%s/bin" % (__dir_cpp_opt, "BUDE-1.2.9", "build")

__dir_getsd_opt = "GetSD"
__dir_getsd_val = "${%s}/%s" % (__dir_cpp_opt, "getSD")

__dir_colsd_opt = "ColourBySD"
__dir_colsd_val = "${%s}/%s" % (__dir_cpp_opt, "colourBySD")

__dir_msrc_opt = "MutationSources"
__dir_mrslt_opt = "mut_rslt"

# Ini file executable options
__exe_bude_opt = "BUDE"
__exe_bude_val = "budeScan"

__exe_bseq_opt = "BudeSeq"
__exe_bseq_val = "bude_sequence"

__exe_colsd_opt = "ColouBySD"
__exe_colsd_val = "colourBySD"

__exe_getsd_opt = "GetSD"
__exe_getsd_val = "getSD"

__exe_fbude_opt = "FullBude"
__exe_fbude_val = "${%s:%s}/${%s}" % (__dir_sect, __dir_bude_opt, __exe_bude_opt)

__exe_fbseq_opt = "FullBseq"
__exe_fbseq_val = "${%s:%s}/${%s}" % (__dir_sect, __dir_bude_opt, __exe_bseq_opt)

__exe_fcolsd_opt = "FullColourSD"
__exe_fcolsd_val = "${%s:%s}/${%s}" % (__dir_sect, __dir_colsd_opt, __exe_colsd_opt)

__exe_fgetsd_opt = "FullGetSD"
__exe_fgetsd_val = "${%s:%s}/${%s}" % (__dir_sect, __dir_getsd_opt, __exe_getsd_opt)

# Ini file names options
__names_b_bmc_opt = "BudeBMC"
__names_b_bmc_val = "budeAlaScan.bemc"
__names_b_ctrl_opt = "BudeCtrl"
__names_b_ctrl_val = "alaScan.bctl"
__names_b_gen_zero_opt = "BudeGenZero"
__names_b_gen_zero_val = "zeroGeneration.bzgn"
__names_b_ff_opt = "BudeForcefield"
__names_b_ff_val = "forcefield.bhff"
__names_b_transf_opt = "BudeTransformations"
__names_b_transf_val = "transformations.bltr"
__names_b_rot_choice_opt = "BudeRotamerChoice"
__names_b_rot_choice_val = "rotamer_choice.brch"
__names_b_rot_lib_opt = "BudeRotamerLibrary"
__names_b_rot_lib_val = "rotamer_library.brot"
__names_ligs_list_opt = "LigandsList"
__names_ligs_list_val = "myLigands.list"
__names_rcpts_list_opt = "ReceptorsList"
__names_rcpts_list_val = "myReceptors.list"
__names_rslt_lig_list_opt = "LigandResultsList"
__names_rslt_lig_list_val = "resultsLigands.list"
__names_rslt_rcpt_list_opt = "ReceptorResultsList"
__names_rslt_rcpt_list_val = "resultsReceptors.list"
__names_lig_sd_opt = "LigandsAvgSD"
__names_lig_sd_val = "ligandsAvgSD.bals"
__names_rcpt_sd_opt = "ReceptorsAvgSD"
__names_rcpt_sd_val = "receptorsAvgSD.bals"

__config_sections = {
    __general_sect:
     {
        __gen_max_auto_num_opt: __gen_max_auto_num_val
     }
    ,
    __exe_sect:
     {
        __exe_bude_opt: __exe_bude_val,
        __exe_bseq_opt: __exe_bseq_val,
        __exe_colsd_opt: __exe_colsd_val,
        __exe_getsd_opt: __exe_getsd_val,
        __exe_fbude_opt: __exe_fbude_val,
        __exe_fbseq_opt: __exe_fbseq_val,
        __exe_fcolsd_opt: __exe_fcolsd_val,
        __exe_fgetsd_opt: __exe_fgetsd_val
     }
    ,
    __dir_sect:
     {
        __dir_wrk_opt: __dir_wrk_val,
        __dir_scan_opt: __dir_scan_val,
        __dir_blibs_opt: __dir_blibs_val,
        __dir_rslt_opt: __dir_rslt_val,
        __dir_src_opt: __dir_src_val,
        __dir_seq_opt: __dir_seq_val,
        __dir_chi_opt: __dir_chi_val,
        __dir_l_src_opt: __dir_l_src_val,
        __dir_l_seq_opt: __dir_l_seq_val,
        __dir_l_blibs_opt: __dir_l_blibs_val,
        __dir_l_rslt_opt: __dir_l_rslt_val,
        __dir_l_chi_opt: __dir_l_chi_val,
        __dir_cpp_opt: __dir_cpp_val,
        __dir_bude_opt: __dir_bude_val,
        __dir_getsd_opt: __dir_getsd_val,
        __dir_colsd_opt: __dir_colsd_val
     }
     ,
    __names_sect:
     {

        __names_b_bmc_opt: __names_b_bmc_val,
        __names_b_ctrl_opt: __names_b_ctrl_val,
        __names_b_gen_zero_opt: __names_b_gen_zero_val,
        __names_b_ff_opt: __names_b_ff_val,
        __names_b_transf_opt: __names_b_transf_val,
        __names_b_rot_choice_opt: __names_b_rot_choice_val,
        __names_b_rot_lib_opt: __names_b_rot_lib_val,
        __names_ligs_list_opt: __names_ligs_list_val,
        __names_rcpts_list_opt: __names_rcpts_list_val,
        __names_rslt_lig_list_opt: __names_rslt_lig_list_val,
        __names_rslt_rcpt_list_opt: __names_rslt_rcpt_list_val,
        __names_lig_sd_opt: __names_lig_sd_val,
        __names_rcpt_sd_opt: __names_rcpt_sd_val
     }

    }


def __create_cfg():
    '''
    Create configuration by adding all sections and options to it.
    '''

    my_cfg = ConfigParser(interpolation=ExtendedInterpolation())
    my_cfg.optionxform = str
    for a_section in __config_sections:
        my_cfg.add_section(a_section)
        for an_option in __config_sections[a_section]:
            my_cfg.set(a_section, an_option, __config_sections[a_section][an_option])

    return my_cfg


def __read_cfg():
    '''
    Read config file, add missing configuration parameters and update
    config file if any section or option was added. 
    '''
    my_cfg = ConfigParser(interpolation=ExtendedInterpolation())
    my_cfg.optionxform = str
    my_cfg.read(__config_file, encoding='utf-8')

    if __did_add_missing_conf(my_cfg):
        copy2(__config_file, (__config_file + ".backup"))
        cfg = open(__config_file, 'w')
        my_cfg.write(cfg)
        cfg.close()

    return my_cfg


def __did_add_missing_conf(my_cfg):
    '''
    Add missing section or option and return true if any change is made.
    '''

    addition_was_made = False

    for a_section in __config_sections:
        if not my_cfg.has_section(a_section):
            my_cfg.add_section(a_section)
            addition_was_made = True
        for an_option in __config_sections[a_section]:
            if not my_cfg.has_option(a_section, an_option):
                my_cfg.set(a_section, an_option, __config_sections[a_section][an_option])
                addition_was_made = True

    return addition_was_made


def create_configuration():
    '''
    Creates configuration or reads it from file.
    '''

    try:
        flags = os.O_EXCL | os.O_WRONLY | os.O_CREAT
        mode = stat.S_IREAD | stat.S_IWRITE | stat.S_IRGRP | stat.S_IROTH
        ini_file = os.open(__config_file, flags, mode)
        os.close(ini_file)

    except OSError as e:
        if e.errno == errno.EEXIST:
            my_config = __read_cfg()
        else:
            raise
    else:
        cfg = open(__config_file, 'w')
        my_config = __create_cfg()
        my_config.write(cfg)
        cfg.close()

    return my_config

