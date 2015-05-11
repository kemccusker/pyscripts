# DDFA.py
import numpy as np
import pandas as pd

basepath='/Users/kelly/Dropbox/projects/Ryan/'
fname = basepath + 'A-B_peak_height.txt'



datdf = pd.read_csv(fname,delimiter="\t")
"""Int64Index: 807 entries, 0 to 806
Data columns (total 9 columns):
Dye/Sample Peak     807  non-null values
Sample File Name    807  non-null values
Marker              305  non-null values
Allele              305  non-null values
Size                787  non-null values
Height              807  non-null values
Area                807  non-null values
Data Point          807  non-null values
Unnamed: 8          0  non-null values
dtypes: float64(2), int64(3), object(4)>
"""
# D/SP          SFN             Marker                          Allele  Size    Height  Area    DataPoint
# "B,65"	A_A03.fsa	_Internal_Marker_Dye_Blue_	91	91.09	2397	13601	1878	


# 1. extract all peak height values within a user-provided range (e.g. size values 100-200):
# @@ question, do you mean use the 'Size' Column to extract the 'Height'?  ANS: YES, have it as user-input
# @@ question, what does it mean when Size has parentheses? ANS: those are alleles. Not Size
# "B,87"	A_A03.fsa	_Internal_Marker_Dye_Blue_	123(1)	122.66	125	452	2058	

# 2. Round size values (floats) up or down to nearest integer value
# 3. For every integer within the user-provided range, check to see if there is a "size" and "height" (if either is missing, both will be). If missing, insert appropriate "size" integer and assign a "height" value of 0.
# @@ question, what do you mean by "appropriate size" here? ANS: within the user specified range, each integer should be accounted for (the Size), so if a size/height is missing, insert the missing index.

# 4. Compare "height" values at each "size" between both samples. If either sample has a 0 "height" at a given "size", set that "height" to 0 in both samples. If easier, can save the values for both samples in one spreadsheet, indicated by different sample file names
# @@ question, isn't that already what you sent me? would prefer A and B samples in separate files actually but only have A right now. ANS: yes already have it
# 5. output both new spreadsheets
# @@ question: what are the two spreadsheets? the 'A' output and the 'diff' output?  ANS: new A, new B, and difference





Aret = datdf.values[datdf['Sample File Name'].values == 'A_A03.fsa']
Bret = datdf.values[datdf['Sample File Name'].values == 'B_B03.fsa']

