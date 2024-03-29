'''
Created on 19 Jun 2017

@author: Amaurys Ávila Ibarra

Transformation and rotation BUDE library.

Syntactically correct but usable for BUDE Docking. 
Not used for Alanine Scanning but mandatory for BUDE.

'''
__version__ = 2.00

trans_rot = '''%BUDE TRANSFORMATIONS LIBRARY
%VERSION 2.00 
# The above 2 lines are required, only the floating point number can vary.
# Any line with a '#' as 1st character is a comment line and ignored, but
# no further comments are allowed once the data entries have started 
# (VERSION 1.00).
# 
# %VERSION 2.00: Will allow comments at any position within the file. Any part
# of a line after a '#' character will considered a comment.It will also
# ignore blank lines.
# 
# Valid transformations are:
# Rotations: "TILT","ROLL","PAN"
# Translations: "XTRANS","YTRANS","ZTRANS"
#
TILT    3
 -10.0
   0.0
  10.0
ROLL    3
 -10.0
   0.0
  10.0
PAN     3 
 -10.0
   0.0
  10.0
XTRANS  3
-1.0
 0.0
 1.0
YTRANS  3
-1.0
 0.0
 1.0
ZTRANS  3
-1.0
 0.0
 1.0
  END
'''

