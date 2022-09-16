'''
Created on 14 Jun 2017

@author: Amaurys Ávila Ibarra
'''

from ntpath import basename
import subprocess, sys

from isambard import ampal


def print_std_error(my_msg):
    print(my_msg, file=sys.stderr)


def run_command(executable, command):
    '''
    Run a command within a shell and pipes standard output and error.
    If return code is not 0 breaks. Otherwise returns the a CompletedProcess object.
    Arguments.
    executable -- The name of the executable that was used in the command.
                  Used only for creating messages.
    command -- Full command as it would have been given in the shell command line.
    
    '''

    executable = basename(executable)

    try:
        proc = subprocess.run([command], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    except Exception as e:
        sys.stderr.write("\nCommand failed: '{}'.\n".format(command))
        sys.stderr.write("Error Info:\n{}\n".format(repr(e)))

    if proc.returncode > 0:
        print("Fatal Error: %s could not run properly. It returns code (%i)" % (executable, proc.returncode), file=sys.stderr)
        if proc.stderr:
            print("%s error message:\n%s\n" % (executable, proc.stderr), file=sys.stderr)
        sys.exit(3)

    return proc


def is_multi_model(my_ampal):
    '''
    Check that we have an Ampal Container (multiple models PDB), if so 
    return true, otherwise return false.
    '''
    if isinstance(my_ampal, ampal.assembly.AmpalContainer):
        return True
    else:
        return False


def create_unique_list(a_list):
    """
    Takes a list and return a list with unique elements.
    
    a_list -- list to remove duplicates from.
    
    
    returns a list of unique elements.
    """
    return list(set(a_list))


def upper_my_list(a_list):
    '''
    Upper-casing all elements of list of strings.
    '''
    return [e.upper() for e in a_list]


def expand_to_single_chars(a_list):
    """
    This function is to expand into single characters values given in the 
    command line like one-letter residue codes or chains ID.    
    
    The elements of the list could be strings or characters.
    
    ie ['DERKH', A] -> ['D', 'E', 'R', 'K', 'H', 'A']
    
    a_list -- List of characters or strings of characters
    
    return a list of characters. 
    """
    expanded_list = []

    for elmnt in a_list:
        if len(elmnt) == 1:
            expanded_list.append(elmnt)
        else:
            for a_char in elmnt:
                expanded_list.append(a_char)

    return expanded_list


def get_legal_chains(my_ampal):
    '''
    Get all chain IDs from the ampal object
    
    my_ampal -- It could be an Ampal Container or a polymer
    '''

    legal_chains = []

    if is_multi_model(my_ampal):
        model = my_ampal[0]
    else:
        model = my_ampal

    for chain in model:
        legal_chains.append(chain.id.upper())

    return legal_chains


def convert_aa(amino, want_longer=False):

    """
    Convert Amino Acids into 3-letter code, 1-letter code or name
    from either the 3-letter code, 1-letter code or name.
    
    want_longer will determine if we want to return the name instead of
    the 3-letter code for a given 1-letter code aa or if we get 
    the 3-letter instead of the 1-letter code for a given name.
    
    ie: convert_aa('A', False) -> 'ALA'
        convert_aa('A', True)  -> 'Alanine'
        
        convert_aa('ALA', False) -> 'A'
        convert_aa('ALA', True)  -> 'Alanine'
        
        convert_aa('Alanine', False) -> 'A'
        convert_aa('Alanine', True)  -> 'ALA'
    
    
    
    amino -- Name or code to convert.
    want_longer -- If we want the longer representation of the amino acid.

    Returns the converted representation of amino, None if amino is not found.
    """

    up_amino = amino.upper()

    amino_acids = {
    "ALA":{"A":"Alanine"},
    "CYS":{"C":"Cysteine"},
    "ASP":{"D":"Aspartic Acid"},
    "GLU":{"E":"Glutamic Acid"},
    "PHE":{"F":"Phenylalanine"},
    "GLY":{"G":"Glycine"},
    "HIS":{"H":"Histidine"},
    "ILE":{"I":"Isoleucine"},
    "LYS":{"K":"Lysine"},
    "LEU":{"L":"Leucine"},
    "MET":{"M":"Methionine"},
    "ASN":{"N":"Asparagine"},
    "PRO":{"P":"Proline"},
    "GLN":{"Q":"Glutamine"},
    "ARG":{"R":"Arginine"},
    "SER":{"S":"Serine"},
    "THR":{"T":"Threonine"},
    "VAL":{"V":"Valine"},
    "TRP":{"W":"Tryptophan"},
    "TYR":{"Y":"Tyrosine"},
    "SEC":{"U":"Selenocysteine"}
    }

    if len(amino) == 1:
        for char3_code in amino_acids:
            char1_code = list(amino_acids[char3_code].keys())[0]

            if char1_code == up_amino and want_longer:
                return amino_acids[char3_code][char1_code]
            elif char1_code == up_amino and not want_longer:
                return  char3_code
        return None
    elif len(amino) == 3:
        if up_amino in amino_acids:
            char1_code = list(amino_acids[up_amino].keys())[0]
            if want_longer:
                return amino_acids[up_amino][char1_code]
            else:
                return char1_code
        else:
            return None
    else:
        for char3_code in amino_acids:
            for char1_code in amino_acids[char3_code]:
                up_name = amino_acids[char3_code][char1_code].upper()
                if up_amino == up_name and want_longer:
                    return char3_code
                elif up_amino == up_name and not want_longer:
                    return char1_code
    return None

