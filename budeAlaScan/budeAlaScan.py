#!/usr/bin/env python3
# encoding: utf-8

'''
budeAlaScan -- Automatically or manually calculate ddG for hot constellations from a PDB.
by doing an Alanine scan on all residues and calculating ddG for each mutant.

It will take the relevant sub-command plus required options for each mode.

Hot constellations will be created from ligand chain(s) in all modes.

All residues in the structure must be standard residues.

Sub commands

auto:
Create all possible hot constellation automatically with all residues
within the cut off distance and ddG greater or equal than the set limit.

manual:
Calculate ddG for the given constellations.

residues:
Create all possible hot constellation with given residues.

scan:
Do alanine scan only. Single mutants, no hot constellations will be created.

@author:     Amaurys Ávila Ibarra

@copyright:  2017 University of Bristol. All rights reserved.

@license:    license

@contact:    Amaurys.AvilaIbarra@bristol.ac.uk
@deffield    updated: Updated
'''

import gc
from ntpath import basename, splitext
import sys, argparse, os
from time import strftime, localtime

from isambard import ampal
from isambard.modelling.scwrl import scwrl_available
from matplotlib.pyplot import show

from budeAlaScan import config
from budeAlaScan.alascan import auto_scan, manual_scan, singles_scan, residues_scan
from budeAlaScan.errorparser.check_input_info import check_input_common_errors_4modes
from budeAlaScan.myutils.check_exe import check_executables
from budeAlaScan.myutils.layout import create_layout
from budeAlaScan.myutils.misce import is_multi_model, expand_to_single_chars
from budeAlaScan.myutils.misce import upper_my_list
from budeAlaScan.plots.plot_results import plot_singles, plot_mutations, write_plot_data

__all__ = []
__version__ = 0.1
__date__ = '2017-04-11'
__updated__ = '2017-04-11'
__sub_commands__ = "auto|manual|residues|scan"
__min_req_python__ = (3, 5)


def main(argv=None):

    try:
        assert sys.version_info >= __min_req_python__
    except Exception as e:
        sys.stderr.write("\nFATAL ERROR: Wrong Python Version:\n\n"
            "We shall NOT continue.\nWe need python {}.{} or greater"
            " for using the budeAlaScan.\n".format(__min_req_python__[0],
                                                   __min_req_python__[1])
            )
        sys.stderr.write("\nWe found this Python Version:\n")
        sys.stderr.write(sys.version)
        sys.stderr.write("\n")
        return 2
#     from budeAlaScan.myutils.test_me import run_my_prog
#     #run_my_prog(7)
#     run_my_prog('1ldb.pdb')
#     return 0

    current_dir = os.getcwd()

    # Update work directory to be current directory
    if current_dir != config.cfg.get(config.dir_sect, config.dir_wrk_opt):
        config.cfg.set(config.dir_sect, config.dir_wrk_opt, current_dir)
#         print("Warning: The work and current directories are not the same.")
#         print("We have updated work directory to (%s)." % (current_dir))
#         print("To remove this warning update the option %s in section '%s'" % (config.dir_wrk_opt, config.dir_sect))
#         print("from ini file (%s)." % (config.config_file))

    if not scwrl_available():
        config.there_is_scwrl = False

    '''Command line options.'''
    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_fulldesc = __import__('__main__').__doc__.split("\n")
    program_longdesc = '\n'.join(program_fulldesc[1:25])

    program_license = '''%s

  Created by Amaurys Ávila Ibarra on %s.
  Copyright 2017 University of Bristol. All rights reserved.

  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.
   
   
 
USAGE
 
%s [%s] [options]

''' % (program_longdesc, str(__date__), program_name, __sub_commands__)

    auto_desc = """
It will do the alanine scanning for all single mutants selecting all residues with 
BUDE's ddG greater or equal than the set limit.

Hot constellations will be created using constellations' size with selected residues if the
distance among their CA is less or equal than the cut off CA-CA distance.
 
It expects the PDB file name, the chain(s) to consider for hot constellations, "ligand chains",
as well as the chain(s) it interacts with, "receptor chain(s)".

The ddG limit, the size of the constellations and the cut off distance from CA to CA can be given
if different from their default values.

Hot constellations will be created only for the ligand chain(s).
"""

    manual_desc = """
It will do the alanine scanning for all single mutants.

It will calculate the ddG for the given hot constellations. A constellation is defined as
a comma separated list of residue numbers preceded by the chain ID. Multiple constellations 
must be space separated.

All residues must be part of the ligand chain(s). 

eg: Two constellations with 3 and 2 residues. First one with residues 3, 4 and 7 chain A
and the second one with residues 21 and 24 from chain B.
[A3,A4,A7 B21,B24]

It expects the PDB file name, the hot constellations to consider, the ligand and receptor chains.
"""
    residues_desc = """
It will do the alanine scanning for all single mutants and will create all possible constellations
of constellation size with the list of given residues.
 
It expects the PDB file name, the residues to create constellations, the ligand and receptor chains.
The size of the constellation can be given if different from default value.

Residues must be given separated by space and should follow the notation chain plus residue number.

The residues must be part of the ligand chain(s).

eg: [A3 A4 A7 B21 B24]
"""

    scan_desc = """
It will do the alanine scanning for all single mutants. No constellation will be created.

It expects the PDB file name, the ligand and receptor chains.

"""

    try:
        # Setup argument parser
        parser = argparse.ArgumentParser(description=program_license, formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument('-V', '--version', action='version', version=program_version_message)
        parser.add_argument('-v', '--verbose', dest="verbose_g", action='store_true', help="Output information during execution [default: %(default)s].")
        parser.add_argument('-t', '--turn-plot-off', dest="no_plot_g", action='store_false', help="Do not show plots [default: %(default)s].")
        parser.add_argument('-i', '--disable-rotamer-correction', dest="no_rot_g", action='store_true', help="No residues will be activated for rotamer correction [default: %(default)s].")

        subparsers = parser.add_subparsers(title="Scan mode sub-commands", description="How to create hot constellations are created.", help="How hot constellation are created.")

        # Auto selection
        parser_auto = subparsers.add_parser('auto', description=auto_desc, help="Create constellations automatically.", formatter_class=argparse.RawDescriptionHelpFormatter)
        parser_auto.add_argument(dest="scan_mode", help=argparse.SUPPRESS, nargs='?', const=True, default='auto')
        parser_auto.add_argument('-v', '--verbose', dest="verbose", action='store_true', help="Output information during execution [default: %(default)s].")
        parser_auto.add_argument('-t', '--turn-plot-off', dest="no_plot", action='store_false', help="Do not show plots [default: %(default)s].")

        # Mandatory auto options.
        mandatory_auto = parser_auto.add_argument_group("Mandatory Arguments")
        mandatory_auto.add_argument("-p", "--pdb-file", dest="pdb_file_name", help="Name of the PDB file.", metavar="myPDB.pdb", required=True)
        mandatory_auto.add_argument("-r", "--receptor-chains", dest="receptor_chains", help="Chain(s) that interact with ligand chain(s).", metavar="A", type=str, nargs='+', required=True)
        mandatory_auto.add_argument("-l", "--ligand-chains", dest="ligand_chains", help="Chain(s) to consider for hot constellations.", metavar="B", type=str, nargs='+', required=True)

        # Default auto options
        default_auto = parser_auto.add_argument_group("Default Arguments")
        default_auto.add_argument("-d", "--deltaDelta-G", dest="deltaDelta_g", help="BUDE ddG selection limit [default: %(default)s].", metavar="5.0", type=float, nargs='?', default=5.0)
        default_auto.add_argument("-z", "--constellation-size", dest="constellation_size", help="Hot constellation size, residues' number [default: %(default)s].", metavar="3", type=int, nargs='?', default=3)
        default_auto.add_argument("-u", "--cut-off", dest="cut_off", const=True, help="Cut off distance (Å) from CA-CA. [default: %(default)s].", metavar="13.0", type=float, nargs='?', default=13.0)
        default_auto.add_argument("-a", "--residues-to-activate", dest="res2activate", help="Residues to activate for rotamer correction. [default: %(default)s].", metavar="D E R K", type=str, nargs='+', default="DERKH")
        default_auto.add_argument('-i', '--disable-rotamer-correction', dest="no_rot_m", action='store_true', help="No residues will be activated for rotamer correction [default: %(default)s].")

        # manual selection
        parser_manual = subparsers.add_parser('manual', description=manual_desc, help="Constellations are given manually, [A1,A3,A6 B4,B6].", formatter_class=argparse.RawDescriptionHelpFormatter)
        parser_manual.add_argument(dest="scan_mode", help=argparse.SUPPRESS, nargs='?', const=True, default='manual')
        parser_manual.add_argument('-v', '--verbose', dest="verbose", action='store_true', help="Output information during execution [default: %(default)s].")
        parser_manual.add_argument('-t', '--turn-plot-off', dest="no_plot", action='store_false', help="Do not show plots [default: %(default)s].")

        # Mandatory manual options.
        mandatory_manual = parser_manual.add_argument_group("Mandatory Arguments")
        mandatory_manual.add_argument("-p", "--pdb-file", dest="pdb_file_name", help="Name of the PDB file.", metavar="myPDB.pdb", required=True)
        mandatory_manual.add_argument("-r", "--receptor-chains", dest="receptor_chains", help="Chain(s) that interact with ligand chain(s).", metavar="A", type=str, nargs='+', required=True)
        mandatory_manual.add_argument("-l", "--ligand-chains", dest="ligand_chains", help="Chain(s) to consider for hot constellations.", metavar="B", nargs='+', type=str, required=True)
        mandatory_manual.add_argument("-c", "--hot-constellations", dest="my_constellations", help="Hot constellations.", metavar="B1,B3,B6", nargs='+', required=True)

        # Default options for full mode.
        default_manual = parser_manual.add_argument_group("Default Arguments")
        default_manual.add_argument("-a", "--residues-to-activate", dest="res2activate", help="Residues to activate for rotamer correction. [default: %(default)s].", metavar="D E", type=str, nargs='+', default="DERKH")
        default_manual.add_argument('-i', '--disable-rotamer-correction', dest="no_rot_m", action='store_true', help="No residues will be activated for rotamer correction [default: %(default)s].")

        # residues selection
        parser_residues = subparsers.add_parser('residues', description=residues_desc, help="Create constellations from residues, [A1 A3 A6 B4 B6].", formatter_class=argparse.RawDescriptionHelpFormatter)
        parser_residues.add_argument(dest="scan_mode", help=argparse.SUPPRESS, nargs='?', const=True, default='residues')
        parser_residues.add_argument('-v', '--verbose', dest="verbose", action='store_true', help="Output information during execution [default: %(default)s].")
        parser_residues.add_argument('-t', '--turn-plot-off', dest="no_plot", action='store_false', help="Do not show plots [default: %(default)s].")

        # Mandatory resiudes options.
        mandatory_residues = parser_residues.add_argument_group("Mandatory Arguments")
        mandatory_residues.add_argument("-p", "--pdb-file", dest="pdb_file_name", help="Name of the PDB file.", metavar="myPDB.pdb", required=True)
        mandatory_residues.add_argument("-r", "--receptor-chains", dest="receptor_chains", help="Chain(s) that interact with ligand chain(s).", metavar="A", type=str, nargs='+', required=True)
        mandatory_residues.add_argument("-l", "--ligand-chains", dest="ligand_chains", help="Chain(s) to consider for hot constellations.", metavar="B", nargs='+', type=str, required=True)
        mandatory_residues.add_argument("-e", "--residues", dest="my_residues", help="Residues to create hot constellations.", metavar="B7", nargs='+', required=True)
        # mandatory_residues.add_argument("-z", "--constellation-size", dest="constellation_size", help="Hot constellation size, number of residues.[default: %(default)s]", metavar="3", type=int, nargs='?', default=3)

        # Default residues options
        default_residues = parser_residues.add_argument_group("Default Arguments")
        default_residues.add_argument("-z", "--constellation-size", dest="constellation_size", help="Hot constellation size, residues' number [default: %(default)s].", metavar="3", type=int, nargs='?', default=3)
        default_residues.add_argument("-a", "--residues-to-activate", dest="res2activate", help="Residues to activate for rotamer correction. [default: %(default)s].", metavar="D E", type=str, nargs='+', default="DERKH")
        default_residues.add_argument('-i', '--disable-rotamer-correction', dest="no_rot_m", action='store_true', help="No residues will be activated for rotamer correction [default: %(default)s].")

        # scan selection
        parser_scan = subparsers.add_parser('scan', description=scan_desc, help="ALA scan single mutants only. No constellations.", formatter_class=argparse.RawDescriptionHelpFormatter)
        parser_scan.add_argument(dest="scan_mode", help=argparse.SUPPRESS, nargs='?', const=True, default='scan')
        parser_scan.add_argument('-v', '--verbose', dest="verbose", action='store_true', help="Output information during execution [default: %(default)s].")
        parser_scan.add_argument('-t', '--turn-plot-off', dest="no_plot", action='store_false', help="Do not show plots [default: %(default)s].")

        # Mandatory manual options.
        mandatory_scan = parser_scan.add_argument_group("Mandatory Arguments")
        mandatory_scan.add_argument("-p", "--pdb-file", dest="pdb_file_name", help="Name of the PDB file.", metavar="myPDB.pdb", required=True)
        mandatory_scan.add_argument("-r", "--receptor-chains", dest="receptor_chains", help="Chain(s) to be consider the receptor docking unit.", metavar="A", type=str, nargs='+', required=True)
        mandatory_scan.add_argument("-l", "--ligand-chains", dest="ligand_chains", help="Chain(s) to be consider the ligand docking unit.", metavar="B", nargs='+', type=str, required=True)

        # Default options for full mode.
        default_scan = parser_scan.add_argument_group("Default Arguments")
        default_scan.add_argument("-a", "--residues-to-activate", dest="res2activate", help="Residues to activate for rotamer correction. [default: %(default)s].", metavar="D E", type=str, nargs='+', default="DERKH")
        default_scan.add_argument('-i', '--disable-rotamer-correction', dest="no_rot_m", action='store_true', help="No residues will be activated for rotamer correction [default: %(default)s].")

        args = parser.parse_args()

    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception as e:
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help\n\n")
        return 2

    if args.verbose or args.verbose_g:
        config.verbose = True

    if not args.no_plot or not args.no_plot_g:
        config.showplots = False

    if args.no_rot_g or args.no_rot_m:
            config.do_rotamer_correction = False

    if config.verbose:
        print("budeAlaScan started on: ", strftime(config.date_fmt, localtime()))
    if not check_executables():
        print("FATAL ERROR: Executable(s) not found.", file=sys.stderr)
        return 2

    if not os.path.isfile(args.pdb_file_name):
        print("FATAL ERROR: The PDB file [%s] cannot be found." % (args.pdb_file_name))
        return 2

    create_layout()

    unpacked_amp = ampal.load_pdb(args.pdb_file_name)
    is_multimodel = is_multi_model(unpacked_amp)

    if is_multimodel:
        model_count = len(unpacked_amp)
        model_pad_size = len(str(model_count)) + 1
    else:
        model_count = 1
        model_pad_size = 0

    # The prefix for the receptor MUST start with 'R' and with 'L' for the ligand.
    rctr_prefix = "Rec"
    lig_prefix = "Lig"

    my_basename = splitext(basename(args.pdb_file_name))[0]

    ligand_chains = sorted(upper_my_list(expand_to_single_chars(args.ligand_chains)))
    receptor_chains = sorted(upper_my_list(expand_to_single_chars(args.receptor_chains)))

    check_input_common_errors_4modes(args.pdb_file_name, receptor_chains,
                                     ligand_chains, unpacked_amp, args.res2activate)

    if args.scan_mode == 'auto':
        if config.verbose:
            print("We are doing an auto hot constellations for PDB: '%s' which has %i model(s)" % (args.pdb_file_name, model_count))
            print("Receptor chain(s): %s" % (args.receptor_chains))
            print("Ligand chain(s): %s" % (args.ligand_chains))
            print("We will use all residues with ddG greater than %.2f kJ/mol." % (args.deltaDelta_g))
            print("CA-CA distance within %.2f %s" % (args.cut_off, config.angstrom_symb))
            print("to create constellations with %i residues." % (args.constellation_size))
        lig_ddgs, rec_ddgs, constellations_delta_g, repack_lig_ddgs = auto_scan(unpacked_amp, is_multimodel, my_basename, rctr_prefix, lig_prefix, args, model_pad_size, receptor_chains, ligand_chains)

    elif args.scan_mode == 'manual':
        if config.verbose:
            print("We are doing a manual hot constellations for PDB: '%s' which has %i model(s)" % (args.pdb_file_name, model_count))
            print("Receptor chain(s): %s" % (args.receptor_chains))
            print("Ligand chain(s): %s" % (args.ligand_chains))
            print("We will use the constellations below:")
            print(args.my_constellations)
        lig_ddgs, rec_ddgs, constellations_delta_g, repack_lig_ddgs = manual_scan(unpacked_amp, is_multimodel, my_basename, rctr_prefix, lig_prefix, args, model_pad_size, receptor_chains, ligand_chains)

    elif args.scan_mode == 'residues':
        if config.verbose:
            print("We are doing hot constellations from give residues for PDB: '%s' which has %i model(s)" % (args.pdb_file_name, model_count))
            print("Receptor chain(s): %s" % (args.receptor_chains))
            print("Ligand chain(s): %s" % (args.ligand_chains))
            print("We will create all possible constellations of size '{}' from residues [{}].".format(args.constellation_size, ", ".join(args.my_residues)))
        lig_ddgs, rec_ddgs, constellations_delta_g, repack_lig_ddgs = residues_scan(unpacked_amp, is_multimodel, my_basename, rctr_prefix, lig_prefix, args, model_pad_size, receptor_chains, ligand_chains)

    elif args.scan_mode == 'scan':
        if config.verbose:
            print("We are doing single alanine scanning for PDB: '%s' which has %i model(s)" % (args.pdb_file_name, model_count))
            print("Receptor chain(s): %s" % (args.receptor_chains))
            print("Ligand chain(s): %s" % (args.ligand_chains))
            print("No hot constellations will be produced.")
        lig_ddgs, rec_ddgs = singles_scan(unpacked_amp, is_multimodel, my_basename, rctr_prefix, lig_prefix, model_pad_size, receptor_chains, ligand_chains)
        constellations_delta_g = None
        repack_lig_ddgs = None

    if config.showplots:
        plot_singles(lig_ddgs, my_basename, model_count)
        wt_ddG = lig_ddgs[0][config.deltaG_key]
        if constellations_delta_g is not None:
            plot_mutations(repack_lig_ddgs, constellations_delta_g,
                           my_basename, model_count, wt_ddG)

    # Writes plot and scan data to files.
    write_plot_data(lig_ddgs, constellations_delta_g, repack_lig_ddgs,
                    my_basename, lig_prefix , model_count, args)
    write_plot_data(rec_ddgs, None, repack_lig_ddgs, my_basename, rctr_prefix,
                    model_count, args)

    if config.verbose:
        print("budeAlaScan successfully finished on: ", strftime(config.date_fmt, localtime()))

    if config.showplots:
        show()

    del lig_ddgs, rec_ddgs, constellations_delta_g
    # END of main
    gc.collect()
    return 0


if __name__ == "__main__":
    # If no arguments is given show help.
    if len(sys.argv) is 1:
        sys.argv.append("-h")

    sys.exit(main())

