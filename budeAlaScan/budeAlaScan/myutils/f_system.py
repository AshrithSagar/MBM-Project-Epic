'''
Created on 5 Jul 2017

@author: Amaurys Ávila Ibarra
'''
import  glob, json
from os import path, makedirs


def create_nested_dir(nested_dir):
    '''
    Create nested directories if they do not exist.
    '''

    if not path.exists(nested_dir):
        makedirs(nested_dir)


def read_file(fname):
    """
    Reads file fname and return its content in list one line per element.

    fname -- file name to read.

    Returns list of lines.
    """
    content = []
    with open(fname, 'r', encoding='utf-8') as itf_file:
        for line in itf_file:
            content.append(line.rstrip('\n\r'))

    return content


def write_file(fname, fcontent):
    """
    Write content to a file. The content could be a single string or a list of
    strings.
    
    fname -- file name to write to.
    fcontent -- file content.
    """
    outf = open(fname, 'w')

    if isinstance(fcontent, list):
        for line in fcontent:
            if isinstance(line, str) and line.endswith('\n'):
                outf.write("%s" % (line))
            else:
                outf.write("%s\n" % (line))
    else:
        outf.write("%s\n" % (fcontent))

    outf.close()

    return


def dump_to_json(data, fname):
    """
    Write data to file fname.

    data -- data to dump to a json file.
    fname -- name of the file to write to.
    """
    f_json = open(fname, "w")
    json.dump(data, f_json, indent=2)
    f_json.close()
    return


def read_json(fname):
    """
    Reads a json file and returns the data:

    fname -- json filename to read

    Returns the data from json file.
    """
    f_json = open(fname, "r")
    my_data = json.load(f_json)
    f_json.close()

    return my_data


def read_csv_file(fname, delim=','):
    """
    Read a comma separated file and return a list of lists. Where each element
    in the parent list is a list of all elements contains in the line.
    
    fname -- file name to read.
    delim -- delimiter to use, comma by default.

    Returns list of list.
    """
    content = []
    with open(fname, 'r', encoding='utf-8') as itf_file:
        for line in itf_file:
            content.append(line.rstrip('\n\r').split(delim))

    return content


def write_csv_file(fname, my_data, delim=','):
    """
    Write my data as comma separated file. It expects list of lists. Where each
    nested list will be a raw of data and each element in the nested list will
    be a column.
    
    fname -- file name to write to.
    my_data -- List of list with data 
    delim -- delimiter to use, comma by default.

    """
    outf = open(fname, 'w')

    for a_row in my_data:
        outf.write("{}\n".format(delim.join(str(col) for col in a_row)))
    outf.close()

    return


def create_backup_fname(fname, dname=''):
    """
    This function will get a file name and count how many are in those adding
    the number at the end of the file.
    
    fname -- File name to backup
    dname -- Directory to where to look for fname.
    
    """
    backup_fname = None

    if dname:
        my_pattern = dname + "/" + fname + "*"
    else:
        my_pattern = fname + "*"

    nmber = len(glob.glob(my_pattern))

    if nmber > 0:
        backup_fname = fname + "_backup_{}".format(nmber)

    return backup_fname

# ##

