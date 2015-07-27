# DDFA.py
import numpy as np
import pandas as pd

basepath='/Users/kelly/Dropbox/projects/Ryan/'
fnameb = basepath + 'B_peak_height.txt'
fnamea = basepath + 'A_peak_height.txt'


sizerange=[100,200]
sizeidx=np.arange(sizerange[0],sizerange[1]+1) # @@@@ NEW question: range inclusive of endpoints?

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


def fillarray2(size, height,sizerange, verb=True):

    """ take size and height arrays, round Size. 
        If an index in sizerange is missing in size, enter 0 for height
    """
    pass
    



#Aret = datdf.values[datdf['Sample File Name'].values == 'A_A03.fsa']
#Bret = datdf.values[datdf['Sample File Name'].values == 'B_B03.fsa']

# A DATASET: ==================
acols=datadf.keys()

#adat = pd.DataFrame(datadf.values[np.logical_and(datadf['Size'] >= sizerange[0], 
#                    datadf['Size'] <= sizerange[1])], columns=acols)
adat = datadf
# Round size data
adat['Size']=datadf['Size'].round()
# Select data based on user Size range:
asel = pd.DataFrame(adat.values[np.logical_and(adat['Size'] >= sizerange[0], 
                                               adat['Size'] <= sizerange[1])], 
                    columns=acols)

# Convert Size to integer vals and use as an index
aidx=asel['Size'].values.astype(int) 
numvals=len(aidx)
# remove duplicate sizes (indices) to check the number of dupes
aidxnodup,reconidx = np.unique(aidx,return_inverse=True) # reconidx reconstructs original
numdups=len(aidx)-len(aidxnodup)

# Get Heights that corresponds to Sizes
heights=asel['Height'].values
# Initialize final height array and final size array, 
#     including room for duplicate size indices
heightary=np.zeros(len(sizeidx)+numdups) #.astype(type(heights[0]))
sizeary=np.arange(sizerange[0],sizerange[1]+1+numdups)

# subtract the first size index to make indices start at zero
#ahidxno=aidxnodup-aidxnodup[0]
ahidx=aidx-aidx[0]

# @@@@@ I don't think this works yet
#   
# Set height data into final height array: 
#     There will be zeros where there are no Size indices
#     Duplicate size indices will have the associated height at that index
heightary[ahidx] =heights.astype(type(heights[0])) # hack to get rid of TypeError

# Set size data into final size array:
#     There should be duplicate size indices in place... @@
sizeary[ahidx] = aidx



# B DATASET: =================
bcols=datbdf.keys()
#bdat = pd.DataFrame(datbdf.values[np.logical_and(datbdf['Size'] >= sizerange[0], datbdf['Size'] <= sizerange[1])], columns=bcols)

bdat = datbdf
bdat['Size']=datbdf['Size'].round()
bsel = pd.DataFrame(bdat.values[np.logical_and(bdat['Size'] >= sizerange[0], bdat['Size'] <= sizerange[1])], columns=bcols)




"""plt.figure()
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
plt.ylabel('Height')"""

# @@@ question: round Size first and then select based on range? Or other way around...
#               ANS: round first, then select

# @@ the round() call on mac is giving an error (appears to be a bug)



def myround(vals):

    rndvals=np.zeros(len(vals))
    for ii in np.arange(0,len(vals)):
        rndvals[ii] = np.round(vals[ii])
        #print vals[ii],rndvals[ii]

    return rndvals



def fillarray(size, height,sizerange, verb=True):

    """ take size and height arrays, round Size. 
        If an index in sizerange is missing in size, enter 0 for height
    """

    sizeidx=np.arange(sizerange[0],sizerange[1]+1) # add extra to get through loop last time
    rndsize=size # already rounded: myround(size)

    if verb:
        print '@@@@@@ rndsize ' + str(rndsize)


    finsize=np.zeros(sizeidx[:-1].shape)
    finheight=np.zeros(sizeidx[:-1].shape)

    datidx=0 # keep track of the index of size data
    ii=0 # keep track of index into sizerange
    for index in sizeidx:
        if verb:
            print datidx,ii,index

        # first check if datidx is past the length of the size array
        if datidx==len(rndsize):
            # We are done. Rest of size and height array is zero
            if verb:
                print '   SIZE ARRAY IS DONE. SET the rest to zero...: ' + str(datidx) 
            break
        # now check if ii is pas the length of the size range array
        if ii==len(sizeidx[:-1]):
            if verb:
                print '   SIZE RANGE is done.... what to do with rest of size array? @@ rndsize.shape ' + str(rndsize.shape) +\
                    ', datidx ' + str(datidx) + ', rest of rndsize array: ' + str(rndsize[datidx:])
            break

        # for each index in size range, check if it exists in size data
        if sizeidx[ii] == rndsize[datidx]:
            if verb:
                print '   Data good: ' + str(datidx) + ', ' + str(ii) + ', ' + str(index) +\
                    ': ' + str(sizeidx[ii]) + ' == ' + str(rndsize[datidx])
            # we are good, the data exists and is good. move it to final array
            finsize[ii] = rndsize[datidx]
            finheight[ii] = height[datidx]
            datidx+=1
            ii+=1
        else:
            # have to loop through rndsize until get to next match, entering zero until then.
            keepgoing=True
            while keepgoing:
                if verb:
                    print '   Entering zero: ' + str(datidx) + ', ' + str(ii) + ', ' + str(index) +\
                        ': ' + str(sizeidx[ii]) + ' != ' + str(rndsize[datidx])

                if  sizeidx[ii] > rndsize[datidx]:
                    print '   We have a repeat Size value! Size:' + str(sizeidx[ii-1]) + ', Height choices: ' +\
                        str(height[ii-1]) + '* and ' + str(height[ii]) + ' (* is saved in return array)'
                    datidx+=1 # move forward in Size array
                    
                else:                                                               
                    finsize[ii] = sizeidx[ii] # put correct size in and move on
                    finheight[ii] = 0
                    #datidx+=1
                    ii+=1

                # check the next index: now does Size match?
                if sizeidx[ii] == rndsize[datidx]:
                    keepgoing=False
                    if verb:
                        print '  YES ' + str(sizeidx[ii]) + ' == ' + str(rndsize[datidx])
                    # jump out of while loop by setting keepgoing to False

    return finsize,finheight





print '======= FILL ARRAYS ================'
#bsizernd=myround(bdat.Size.values)

#bsizeadj,bheightadj=fillarray(bdat.Size.values, bdat.Height, sizerange)

#asizernd=myround(adat.Size.values)
#asizeadj,aheightadj=fillarray(adat.Size.values, adat.Height, sizerange)


asizeadj,aheightadj=fillarray(asel.Size.values, asel.Height, sizerange)

bsizeadj,bheightadj=fillarray(bsel.Size.values, bsel.Height, sizerange)



# ==================== FIGURES ========================

fig,axs=plt.subplots(3,1)
fig.set_size_inches(8,10)
ax=axs[0]
ax.plot(bsel.Size,bsel.Height,marker='o',color='r')
ax.plot(asel.Size,asel.Height,marker='.',color='b')
#ax.set_xlabel('Size')
ax.set_ylabel('Height')
ax.set_title('Size Range (' + str(sizerange[0]) + '-' + str(sizerange[1]) + ')')
ax.set_xlim(sizerange)
ax.legend(('B','A'))

ax=axs[1]
ax.plot(bsizeadj,bheightadj,marker='o',color='r')
ax.plot(asizeadj,aheightadj,marker='.',color='b')
#ax.set_xlabel('Size')
ax.set_ylabel('Adjusted Height')
ax.set_xlim(sizerange)
ax.legend(('B','A'))


ax=axs[2]
# @@@ won't work until adjust Size vals with zeros.
ax.plot(asizeadj,aheightadj-bheightadj,marker='.',color='k')
ax.axhline(y=0,color='k')
ax.set_xlabel('Size')
ax.set_ylabel('Height')
ax.set_title('Difference (A-B)')
ax.set_xlim(sizerange)
fig.savefig('peakal.pdf')

