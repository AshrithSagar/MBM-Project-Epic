#!/usr/bin/env python3
# encoding: utf-8

'''
Created on 21 May 2019

@author: Amaurys Ávila Ibarra

BUDE rotamer choice library. 

Optional library for BUDE that it is required for mutagenesis.

'''
__version__ = 2.00

rotamer_choice = """%BUDE ROTAMER CHOICE
%VERSION 2.00
# Remember that the 2 previous lines are required,
# This file is dependant on the rotamer library.
# This file allows you to specify which rotamers of a particular amino acid 
# you want to use for the energy calculations.
#
# %VERSION 2.00: The first line of data must have two intergers. The number of
# amino acids and the number of rotamers.The rotamers must not exceed the 
# MAX_ROTAMER_CHOICE_NUMBER, currently set to 200.
# Amin oacids must not have 1 in any rotamer higher than those declare in the
# rotamer library.
#
# 1 means use the rotamer, 0 means do not use the rotamer
#
# %VERSION 2.00: Will allow comments at any position within the file. Any part
# of a line after a '#' character will considered a comment.It will also
# ignore blank lines.
#
# %VERSION 2.00: The tag END will denote the end of the library. No further 
# lines will be considered. The line containing the END tag can only have the 
# tag and blank characters or a comment.
#
# If the END tag is omitted the file will be process till the EOF and a warning 
# will be issued.
# 
#    1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20    
20 20
ALA  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    
ARG  1  1  0  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0      
ASN  1  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0     
ASP  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0     
CYS  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0     
GLN  1  1  1  1  1  1  1  1  1  1  1  1  1  0  0  0  0  0  0  0     
GLU  1  1  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0    
HIS  1  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0    
ILE  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    
LEU  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    
LYS  1  1  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0    
MET  1  1  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0    
PHE  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    
PRO  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0     
SER  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    
THR  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    
TRP  1  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0    
TYR  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    
VAL  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    
GLY  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    
END  
"""

