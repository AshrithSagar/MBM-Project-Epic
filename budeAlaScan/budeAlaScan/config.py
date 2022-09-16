'''

Configuration module.

Created on 20 Jun 2017

@author: Amaurys Ávila Ibarra
'''

from . import initialize


cfg = initialize.create_configuration()

config_file = initialize.__config_file
exe_sect = initialize.__exe_sect
dir_sect = initialize.__dir_sect
names_sect = initialize.__names_sect
general_sect = initialize.__general_sect

angstrom_symb = u'\u212B'

# Indices from BUDE results in the order they appear blank space separated
# after splitting a line.
# Index Number Name Chain     InterDG    InterDDG  NormTerDDG     IntraDG    IntraDDG  NormTraDDG ChainAtoms     SD
#     1    340  ALA     B   -326.3490      0.0000      0.0000  -1264.1267      0.0000      0.0000          0     0.0000
#
# b_idx -- Index from BUDE, unique incremental integer.
# b_rNum_idx -- Index for the residue number.
# b_rNam_idx -- Index for the residue name.
# b_ch_idx -- Index for the chain ID.
# b_rddG_idx -- Index for the ddG value.
# b_sch_idx -- Index for the sidechain size.
# b_sd_idx -- Index for the standard deviation.

b_idx = 0
b_rNum_idx = 1
b_rNam_idx = 2
b_ch_idx = 3
b_rddG_idx = 5
b_sch_idx = 10
b_sd_idx = 11

# Indices for the tuples of the results we process.
r_rNum_idx = 0
r_rNam_idx = 1
r_ch_idx = 2
r_rddG_idx = 3
r_sch_idx = 4
r_sd_idx = 5

# Key for deltaG in dictionary
deltaG_key = 'dG'

# Re-packing suffix
rpck_suffix = "repacked"

# Replot directory
replot_dir = "replot"

# Directory for file that exceeds Max number constellations
extra_const_dir = "extra_const"

# If there is scwrl
there_is_scwrl = True

# Flag to output info to the standard output.
verbose = False
# Flag to show plots
showplots = True
# Time format for verbose output.
time_fmt = "%H:%M:%S"
# date format for verbose output
date_fmt = "%c"

gen_max_auto_num_opt = initialize.__gen_max_auto_num_opt
dir_wrk_opt = initialize.__dir_wrk_opt
dir_scan_opt = initialize.__dir_scan_opt
dir_blibs_opt = initialize.__dir_blibs_opt
dir_rslt_opt = initialize.__dir_rslt_opt
dir_src_opt = initialize.__dir_src_opt
dir_chi_opt = initialize.__dir_chi_opt
dir_l_src_opt = initialize.__dir_l_src_opt
dir_l_seq_opt = initialize.__dir_l_seq_opt
dir_l_blibs_opt = initialize.__dir_l_blibs_opt
dir_l_rslt_opt = initialize.__dir_l_rslt_opt
dir_l_chi_opt = initialize.__dir_l_chi_opt
dir_cpp_opt = initialize.__dir_cpp_opt
dir_bude_opt = initialize.__dir_bude_opt
dir_seq_opt = initialize.__dir_seq_opt
dir_getsd_opt = initialize.__dir_getsd_opt
dir_colsd_opt = initialize.__dir_colsd_opt
dir_msrc_opt = initialize.__dir_msrc_opt
dir_mrslt_opt = initialize.__dir_mrslt_opt
exe_bude_opt = initialize.__exe_bude_opt
exe_bseq_opt = initialize.__exe_bseq_opt
exe_colsd_opt = initialize.__exe_colsd_opt
exe_getsd_opt = initialize.__exe_getsd_opt
exe_fbude_opt = initialize.__exe_fbude_opt
exe_fbseq_opt = initialize.__exe_fbseq_opt
exe_fcolsd_opt = initialize.__exe_fcolsd_opt
exe_fgetsd_opt = initialize.__exe_fgetsd_opt
names_b_bmc_opt = initialize.__names_b_bmc_opt
names_b_ctrl_opt = initialize.__names_b_ctrl_opt
names_b_gen_zero_opt = initialize.__names_b_gen_zero_opt
names_b_ff_opt = initialize.__names_b_ff_opt
names_b_transf_opt = initialize.__names_b_transf_opt
names_b_rchoice_opt = initialize.__names_b_rot_choice_opt
names_b_rlib_opt = initialize.__names_b_rot_lib_opt
names_ligs_list_opt = initialize.__names_ligs_list_opt
names_rcpts_list_opt = initialize.__names_rcpts_list_opt
names_rslt_lig_list_opt = initialize.__names_rslt_lig_list_opt
names_rslt_rcpt_list_opt = initialize.__names_rslt_rcpt_list_opt
names_lig_sd_opt = initialize.__names_lig_sd_opt
names_rcpt_sd_opt = initialize.__names_rcpt_sd_opt

# Whether we want to do residue rotmamer correction or not.
do_rotamer_correction = True

# Residues to activate to add rotamers for correction. Default charge resiudes.
res2activate = ("D", "E", "H", "K", "R")

# Legal one letter code for the 20 protein amino acids.
legal_aa = ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P",
            "Q", "R", "S", "T", "V", "W", "Y")

rec_bseq_fname = None
lig_bseq_fname = None

