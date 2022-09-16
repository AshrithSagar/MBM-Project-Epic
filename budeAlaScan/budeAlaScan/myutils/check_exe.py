'''
Created on 21 Jun 2017

@author: Amaurys Ávila Ibarra


Checking that all required executables are in path or we have absolute path for them.

'''

from os import path
import shutil

from budeAlaScan import config


def found_cmd(cmd):
    '''
    Check if a command exist, if so return True
    otherwise return false.
    '''
    my_exe = shutil.which(cmd)

    if my_exe is None:
        return False
    else:
        return True


def find_exe_locations(cmd, full_cmd, is_codedir_empty, exe_key):
    '''
    Check that the executable is the path or if relative update to absolute path.
    If the executable is found in either location return True
    otherwise return False.
    '''

    if found_cmd(cmd):
        # Command is in the PATH
        return True
    elif found_cmd(full_cmd):
        # If full_cmd is relative let's change it into full path.
        config.cfg.set(config.exe_sect, exe_key, path.abspath(full_cmd))
        return True

    if not found_cmd(cmd) and not found_cmd(full_cmd):
        if is_codedir_empty:
            # Is it a local directory with no CPP code directory? Let's remove '/'
            local_cmd = full_cmd[1:]
            if found_cmd(local_cmd):
                config.cfg.set(config.exe_sect, exe_key, path.abspath(local_cmd))
                return True
            else:
                print("Error:\nCommand [%s] was not found in the PATH." % (cmd))
                print("Neither we could find [%s] locally.\n" % (local_cmd))
                return False

        print("Error:\nCommand [%s] was not found in the PATH." % (cmd))
        print("Neither we could find [%s] locally.\n" % (full_cmd))
        return False


def check_executables():
    '''
    Return True if all executables are present
    otherwise return False.
    '''
    exe_ok = True

    code_dir = config.cfg.get(config.dir_sect, config.dir_cpp_opt)

    bude_cmd = config.cfg.get(config.exe_sect, config.exe_bude_opt)
    get_sd_cmd = config.cfg.get(config.exe_sect, config.exe_getsd_opt)
    colour_cmd = config.cfg.get(config.exe_sect, config.exe_colsd_opt)
    bseq_cmd = config.cfg.get(config.exe_sect, config.exe_bseq_opt)

    bude_full = config.cfg.get(config.exe_sect, config.exe_fbude_opt)
    bseq_full = config.cfg.get(config.exe_sect, config.exe_fbseq_opt)

    get_sd_full = config.cfg.get(config.exe_sect, config.exe_fgetsd_opt)
    colour_full = config.cfg.get(config.exe_sect, config.exe_fcolsd_opt)

    # Assign False to show all missing executables as oppose to return.
    if not find_exe_locations(bude_cmd, bude_full, not code_dir, config.exe_bude_opt):
        exe_ok = False

    if not find_exe_locations(bseq_cmd, bseq_full, not code_dir, config.exe_bseq_opt):
        exe_ok = False

    if not find_exe_locations(get_sd_cmd, get_sd_full, not code_dir, config.exe_getsd_opt):
        exe_ok = False

    if not find_exe_locations(colour_cmd, colour_full, not code_dir, config.exe_colsd_opt):
        exe_ok = False

    return exe_ok

