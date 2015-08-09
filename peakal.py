""" peakal.py
	Align peaks of samples in specified input file in delimited text format

	Aug 7, 2015
"""

import DDFA as ddfa
import sys

printtofile=False
plotfig = True
verb = True


# Set up inputs
#basepath='./Ryan/'
#ext = '.txt'
#fname = basepath + 'B_peak_height'
#fnamea = basepath + 'A_peak_height'
#fname1 = basepath + 'set1'
#fname2 = basepath + 'set2'
#fname3 = basepath + 'set3'
#fnameh = basepath + 'hnsdata'
#fnamebug = basepath + 'hns_c_03_test'

#delimiter='\t'

#fin = fname+ext
#fout = fname + '_out' + ext

#sizelims=100,200

# initialize
delimiter = '\t'


# Command line args:

if len(sys.argv)>1:
    fin = sys.argv[1]
    sizelims = sys.argv[2]
    if len(sys.argv)==4:
        fout = sys.argv[3]
    else:
        fout = None
else:
    # @@@@ check for python 3, then the func is just input()
    fin = raw_input('Enter input filename (must be in current dir, or provide full path): ')
    delimiter = raw_input('Input filename delimiter (default ''\t''): ' )
    sizelims = raw_input('Enter size range (end point inclusive). E.g. 100,200: ')
    fout = raw_input('Enter output filename. Press enter for default: ')


if fout=='': fout = None
if delimiter=='': delimiter = '\t'



ddfa.test_main(fin, sizelims, outfile=fout, plotfig=plotfig, verb=verb, printtofile=printtofile)

# in HNS file:
# check C_C03.fsa
# sizes ~ 150-156
# for some reason, duplicates found is incorrect:
#   Duplicates found: [122 156 200 150 197 200]
