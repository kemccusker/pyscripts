# DDFA.py
import numpy as np
import pandas as pd

basepath='/Users/kelly/Dropbox/projects/Ryan/'
fnameb = basepath + 'B_peak_height.txt'
fnamea = basepath + 'A_peak_height.txt'


sizerange=[100,200]

datadf = pd.read_csv(fnamea,delimiter="\t")
datbdf = pd.read_csv(fnameb,delimiter="\t")
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





#Aret = datdf.values[datdf['Sample File Name'].values == 'A_A03.fsa']
#Bret = datdf.values[datdf['Sample File Name'].values == 'B_B03.fsa']


acols=datadf.keys()
bcols=datbdf.keys()

adat = pd.DataFrame(datadf.values[np.logical_and(datadf['Size'] >= sizerange[0], datadf['Size'] <= sizerange[1])], columns=acols)
bdat = pd.DataFrame(datbdf.values[np.logical_and(datbdf['Size'] >= sizerange[0], datbdf['Size'] <= sizerange[1])], columns=bcols)


plt.figure()
plt.plot(bdat.Size,marker='o',color='r')
plt.plot(adat.Size,marker='.',color='b')
plt.title('Size')


plt.figure() 
plt.plot(bdat.Height,marker='o',color='r')
plt.plot(adat.Height,marker='.',color='b')
plt.title('Height')

plt.figure()
plt.plot(bdat.Size,bdat.Height,marker='o',color='r')
plt.plot(adat.Size,adat.Height,marker='.',color='b')
plt.xlabel('Size')
plt.ylabel('Height')

# @@@ question: round Size first and then select based on range? Or other way around...

# @@ the round() call on mac is giving an error (appears to be a bug)



def myround(vals):

    rndvals=np.zeros(len(vals))
    for ii in np.arange(0,len(vals)):
        rndvals[ii] = np.round(vals[ii])
        print vals[ii],rndvals[ii]

    return rndvals



def fillarray(size, height,sizerange):

    """ take size and height arrays, round Size. 
        If an index in sizerange is missing in size, enter 0 for height
    """

    sizeidx=np.arange(sizerange[0],sizerange[1]+1)
    rndsize=myround(size)
    
    datidx=0 # keep track of the index of size data
    for ii,index in enumerate(sizeidx):
        # for each index in size range, check if it exists in size data
        if sizeidx[ii] == rndsize[datidx]:
            # we are good, the data exists and is good. move it to final array
            finsize[ii] = rndsize[datidx]
        else:
            # have to loop through rndsize until get to next match, entering zero until then.
            pass

    return sizeadj,heightadj


fig,axs=plt.subplots(2,1)
ax=axs[0]
ax.plot(myround(bdat.Size.values),bdat.Height,marker='o',color='r')
ax.plot(myround(adat.Size.values),adat.Height,marker='.',color='b')
ax.set_xlabel('Size')
ax.set_ylabel('Height')
ax.set_title('Rounded Size')
ax.legend(('B','A'))

ax=axs[1]
# @@@ won't work until adjust Size vals with zeros.
#ax.plot(myround(adat.Size.values),adat.Height-bdat.Height,marker='.',color='k')
ax.set_xlabel('Size')
ax.set_ylabel('Height')
ax.set_title('Difference')
