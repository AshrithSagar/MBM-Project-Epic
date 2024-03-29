'''
Created on 19 Jun 2017

@author: Amaurys Ávila Ibarra

This will contain the BUDE EMC generation library.

Syntactically correct but usable for BUDE Docking. 
Not used for Alanine Scanning but mandatory for BUDE.
'''
__version__ = 3.00

bemc = '''%BUDE EMC FILE
%VERSION  3.00
# Remember that the 2 previous lines are required,
# any other lines starting with # are comments.
# The number of generations is defined by an integer in the first
# 20 characters of the next, non-comment, line.
# No more comments allowed once data has started.
#
# %VERSION 2.00: Will allow comments at any position within the file. Any part
# of a line after a '#' character will considered a comment.It will also
# ignore blank lines.
#
#____________________________________________________________________________________________________________________________________________________
# int    | int     | int     | int     | int      | int      | int      | char     | int     | double   | double   | short int| Boolean  | Boolean |
# use    | select  | print   | print   | create   | print    | print    | choose   | choose  | reserved | reserved | Do       | Select   | Print   |
# Seed   | Parents | Parents | Parents | Children | Children | Children | Mutation | Max Mut | param2   | param3   | Mini GA  | Unique   | Unique  |
#        | Number  | Energy  | Coord   | Number   | Energy   | Coord    | Method   | States  |          |          |          | Parent   | Children|
#________|_________|________ |_________|__________|__________|__________|__________|_________|__________|__________|__________|__________|_________|_
 1   
#____________________________________________________________________________________________________________________________________________________
26182060       100        5          0       5000        10          10       2          4        2.00      3.00        0        true       false 
END

'''
