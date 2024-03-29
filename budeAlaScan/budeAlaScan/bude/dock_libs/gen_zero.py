'''
Created on 19 Jun 2017

@author: Amaurys Ávila Ibarra

Generation Zero for BUDE Docking

Syntactically correct but usable for BUDE Docking. 
Not used for Alanine Scanning but mandatory for BUDE.
'''
__version__ = 2.00

gen_zero = '''%BUDE TRANSFORMATION DESCRIPTOR FILE
%VERSION  2.0
#
# Rotations: "TILT","ROLL","PAN"
# Translations: "XTRANS","YTRANS","ZTRANS"
#
# The order of the states should not be altered.
# it MUST always be
# TILT ROLL PAN XTRANS YTRANS ZTRANS
#
#
#     val1      val2      val3
#        v         v         v
#
             7             6            10
      1      1      1      1      2      2
      3      2      1      1      3      1
      2      1      1      3      3      1
      3      1      1      2      2      1
      1      2      1      1      3      2
      2      3      3      1      1      1
      3      1      2      2      1      1
      2      2      2      3      1      1
      1      2      2      2      1      1
      2      2      2      1      2      1
END

'''
