'''

Repack all chains for the ampal object accordingly.

Created on 30 Jun 2017

@author: Amaurys Ávila Ibarra
'''
from copy import deepcopy

import ampal
from isambard.modelling import scwrl


def getAssemblyFullSequence(my_assembly):

    full_seq = ""

    for poli_pept in my_assembly:
        full_seq += poli_pept.sequence

    return full_seq

def repackAmpal(my_ampal, is_multimodel, my_sequence=None):
    """
    Function for repacking an ampal object. If we have multiple models we will
    repack one at the time. Each model must have identical sequence.
    
    my_ampal -- Ampal object from PDB, container or assembly.
    my_sequence -- Sequence for the object, list of strings, one per chain.
    
    return repacked ampal object.
    """

    rpck_ma = deepcopy(my_ampal)

    if is_multimodel:
        if not my_sequence:
            my_sequence = my_ampal[0].sequences
        repacked_ampal = ampal.AmpalContainer()
        for a_mdl in rpck_ma:
            tmp_mdl = scwrl.pack_side_chains_scwrl(a_mdl, my_sequence)
            repacked_ampal.append(tmp_mdl)
    else:
        if not my_sequence:
            my_sequence = my_ampal.sequences
        repacked_ampal = scwrl.pack_side_chains_scwrl(rpck_ma, my_sequence)

    return repacked_ampal

#####
