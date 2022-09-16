'''
Created on 19 Jun 2017

@author: Amaurys Ávila Ibarra

BUDE control file common values.

For more information use BUDE command and option '-h'
or '-B' for a full control file template.

'''

from budeAlaScan import config

__all_ = [ "bude_ctrl" ]

libs_dir = config.cfg.get(config.dir_sect, config.dir_blibs_opt)
rslts_dir = config.cfg.get(config.dir_sect, config.dir_rslt_opt)
zero_gen_name = config.cfg.get(config.names_sect, config.names_b_gen_zero_opt)
transf_name = config.cfg.get(config.names_sect, config.names_b_transf_opt)
ff_name = config.cfg.get(config.names_sect, config.names_b_ff_opt)
emc_name = config.cfg.get(config.names_sect, config.names_b_bmc_opt)
rot_choice_name = config.cfg.get(config.names_sect, config.names_b_rchoice_opt)
rot_lib_name = config.cfg.get(config.names_sect, config.names_b_rlib_opt)

bude_ctrl = '''# This is the control file. It should comply with the following format.
# Every line that begins with '#' would be treated as comment.
# In-line comments are allowed. After a '#' all will be considered a comment.
#
# Any empty line or with white spaces will be ignored.
# This file may contain none, any or all options. If an option is repeated
# the value of the last occurrence will overwrite any previous one.
# ie.
# option1 value1.
# option1 value2.
# The value of option1 will be value2.
#
# The END tag will denote the ending of the significant part of the file.
# If the tag is omitted the file will be processed till the EOF

# This flag will set whether we want to run 
# a docking that will use the centre of mass of the molecule
# as the reference to apply the transformations.
# Values (true, false).
# Default value 'true'
# command line option(s): -l --use-mass-centre
--use-mass-centre false

# ----------------8<-------------------8<---------------------
# GENERATION ZERO.
# Zero generation file name.
# must exist if to be used
# If --zero-generation-filename is given then
# --generation-zero-descriptor-number MUST be 0 or empty
# --exhaustive-generation MUST  be set to false
# %%BUDE TRANSFORMATION DESCRIPTOR FILE
# default value ''.
# file extension: .bzgn
# command line option(s): -Z --zero-generation-filename
--zero-generation-filename {0}/{1}
# ----------------8<-------------------8<---------------------
# FILES OPTIONS.
# The three options below will affect the naming
# of the files being output by BUDE
#
# Output file's Base. (file.ext -> file)
# This option is used to generated the output files.
# default value 'compound'.
# command line option(s): -R --output-file-base
--output-file-base {2}

# ----------------8<-------------------8<---------------------
# RECEPTOR FILES.
# File Options for the receptor and Docking Grid Centres,
# only the coordinates file is mandatory.
#
# Receptor coordinates file name.
# Must exist.
# PDB, Mol2 format.
# Or a file with a list of file names. One file name per line
# @<TRIPOS>MOLECULE
# %%BUDE LIST FILE
# file extension: .blst
# All files in the list must be in PDB or Mol2 format
# default value 'receptor.pdb'.
# command line option(s): -P --receptor-coordinates-filename
# --receptor-coordinates-filename # THIS must be given in the command line.

# file extension: .bseq
# command line option(s): -S --receptor-sequence-filename
--receptor-sequence-filename 

# ----------------8<-------------------8<---------------------
# LIGAND FILES.
# File Options for the Ligand, only the coordinates file is mandatory.
#
# Ligand coordinates file name.
# Must exist.
# PDB, Mol2 format.
# Or a file with a list of file names. One file name per line
# @<TRIPOS>MOLECULE
# %%BUDE LIST FILE
# file extension: .blst
# All files in the list must be in PDB or Mol2 format
# default value 'ligand.pdb'.
# command line option(s): -I --ligand-coordinates-filename
# --ligand-coordinates-filename # THIS must be given in the command line.

# file extension: .bseq
# command line option(s): -Q --ligand-sequence-filename
#--ligand-sequence-filename # THIS must be given in the command line.
# ----------------8<-------------------8<---------------------
# MANDATORY BUDE LIBRARIES.
# Mandatory BUDE libraries for all docking cases.

# Ligand's Translations - Rotations file name
# Must exist.
# %%BUDE TRANSFORMATIONS LIBRARY
# default value 'ligand_trans_rot.bltr'.
# command line option(s): -T --transformations-filename
--transformations-filename {0}/{3}

# Forcefield file name.
# Must exist.
# %%BUDE HEAVY ATOM AMINO-ACID FORCEFIELD LIBRARY
# default value 'force_field.bhff'.
# command line option(s): -F --forcefield-filename
--forcefield-filename {0}/{4}

# EMC evolution file name.
# Must exist.
# %%BUDE EMC FILE
# default value 'compound.bemc'.
# command line option(s): -V --emc-filename
--emc-filename {0}/{5}

# ----------------8<-------------------8<---------------------
# Activate Residues Mode. Values 00 - 22
# This is a two digits base 3 value. The first digit is for the
# receptor and the second for the ligand.
# 0 - means do not activate.
# 1 - means activate by sequence file.
# 2 - means activate by surface file.
# ie 12 - means activate receptor by sequence  and the ligand by surface.
# Default value '00'.
# command line option(s): -a --residue-activation-mode
--residue-activation-mode 11

# OPTIONAL BUDE LIBRARIES.
# Optional BUDE libraries, but some of them
# might become mandatory depending on the docking case.

# Rotamer library file name.
# must exist if to be used
# %BUDE ROTAMER LIBRARY
# default value ''.
# file extension: .brot
# command line option(s): -L --rotamer-library-filename
--rotamer-library-filename {0}/{7}

# Rotamer choice file name.
# must exist if to be used
# %BUDE ROTAMER CHOICE
# default value ''.
# file extension: .brch
# command line option(s): -C --rotamer-choice-filename
--rotamer-choice-filename {0}/{6}

# ----------------8<-------------------8<---------------------

END # Any information after this tag will be ignored.

'''.format(libs_dir, zero_gen_name, rslts_dir, transf_name, ff_name, emc_name,
           rot_choice_name, rot_lib_name)

