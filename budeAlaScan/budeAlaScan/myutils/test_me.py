'''

A module for doing different tests.
This could be safely deleted once all is working.

Created on 30 Jun 2017

@author: Amaurys Ávila Ibarra
'''

from isambard import ampal

from budeAlaScan.myutils.misce import is_multi_model
from budeAlaScan.myutils.repack import repackAmpal


def run_my_prog(value):
    '''
    This just a small function for testing small portions of the code
    This does not do any significant part of the work.
    
    value -- It is a value which could change depending on the test required.
    '''

    my_ampal = ampal.convert_pdb_to_ampal(value)
    is_multimodel = is_multi_model(my_ampal)

    repack_ampal = repackAmpal(my_ampal, is_multimodel)

    print("\n\n", type(repack_ampal), "\n\n")

    del(repack_ampal)

    if value == 1:
        var1 = "true"
    else:
        var1 = "false"

    print("Var1: ", var1)
    a_dict = {}
    for i in range(1, 10):
        a_list = []
        for k in range(1, 5):
            a_list.append(k)
        a_dict[i] = a_list

    print("A LIST: ", a_list)

    for i in range(1, 10):
        print(i, ": ", a_dict[i])

    return
