""" DDFA.py


    Data example:

    Int64Index: 807 entries, 0 to 806
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


 D/SP          SFN             Marker                          Allele  Size    Height  Area    DataPoint
 "B,65"	       A_A03.fsa	_Internal_Marker_Dye_Blue_	91	91.09	2397	13601	1878	

    ALGORITHM:
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

"""
import numpy as np
import pandas as pd
import cccmautils as cutl
import copy as copy


def process_sample(datdf, sizelims, selcol='Size', heightcol='Height',
                   delimiter="\t", plotfig=True, verb=True, printtofile=False):
    """ do one sample
    
                sizelims: range of sizes to process. Endpoints inclusive

                returns:  ((sizeary, heightary), 
                           (sizearyunq, heightaryunq),
                           (sizeraw, heightraw),
                           sizerawary)
                     where sizerawary are the raw sizes with NaNs as filler
    """

    sizerange=np.arange(sizelims[0],sizelims[1]+1) 

    #datdf = pd.read_csv(fname,delimiter=delimiter)
    fname = datdf['Sample File Name'].values[0] # should only be one sample file name

    cols = datdf.keys()
    # save original size
    sizerawall = datdf[selcol].astype(np.float) 
    # save raw selection
    sizeraw = sizerawall[np.logical_and(sizerawall.round() >= sizelims[0],
                                        sizerawall.round() <= sizelims[1])].values

    # Round size data
    datdf['Size'] = datdf['Size'].astype(np.float).round()

    # Select data based on user Size (selcol) range:
    selbool=np.array(np.logical_and(datdf[selcol] >= sizelims[0],datdf[selcol] <= sizelims[1]))
    seldf = pd.DataFrame(datdf.values[selbool], columns=cols)

    # Convert Size (selcol) to integer vals and use as an index
    rawidx=seldf[selcol].values.astype(int) 
    numsize=len(rawidx)


    # Find duplicates:
    # Use dupe indices to insert dupe data later - these are relative to size data
    # (this next line works b/c rawidx is sorted, otherwise, sort first)
    # Find elements that are the same as the previous element
    dupes=rawidx[rawidx[1:] == rawidx[:-1]]
    if verb:
        print 'Duplicates found: ' + str(dupes)
        print 'rawidx[1:] == rawidx[:-1] ' + str(rawidx[1:] == rawidx[:-1])

    # remove duplicate sizes (indices) so just have unique sizes/indices
    # reconidx reconstructs original, retidx are the indices that result in unique array
    # Will use retidx to set height data associatd w/ unique size indices
    unqidx, retidx, reconidx = np.unique(rawidx,return_index=True,return_inverse=True) 

    # dupeinsert indices are relative to original raw array inex. 
    # Use to access raw height
    dupesinsert=np.zeros((len(dupes)))
    for dii,dup in enumerate(dupes):
        # get index of duplicate size
        dupesinsert[dii] = np.int(cutl.find_nearest(rawidx,dup))

    # Get Heights that corresponds to Sizes
    heightraw=seldf[heightcol].values


    # Initialize final height array and final size array 
    heightary=np.zeros(len(sizerange)) # height should be zero unless there is a size
    sizeary=np.arange(sizelims[0],sizelims[1]+1) # sizes are indices w/in user range
    # Initialize final raw size array (with NaNs as filler)
    sizerawary = np.ones(len(sizerange))*np.nan
    # Initialize 'notes/comments' array
    notesary = ('',)*len(sizerange)

    # Put raw height data into final height array where unique sizes exist
    heightraw = heightraw.astype(type(heightraw[0])) # make sure types match to avoid TypeError
    print '@@@ type(heightraw) ' + str(type(heightraw))
    heightary[unqidx-sizelims[0]] = heightraw[retidx]

    # Put raw size data into final raw size array for unique (rounded) sizes
    sizerawary[unqidx-sizelims[0]] = sizeraw[retidx]

    # Save arrays WITHOUT DUPLICATES:
    sizearyunq = copy.copy(sizeary)
    heightaryunq = copy.copy(heightary)
    sizerawaryunq = copy.copy(sizerawary)


    # ADD DUPLICATES BACK IN: 
    #  convert to list to do insertion
    heightl = list(heightary)
    sizel = list(sizeary)
    sizerawl = list(sizerawary)
    notesl = list(notesary)
    dii=0
    incr=0 # keep track of how many elements we add.
    for dupii in dupes: # dupe size values (also indices into full size array)

        # list.insert(index,value)
        # inserts the value before the index (so it becomes that index). 
        #   dupii is index of first duplicate,
        #   add one to get to next spot, plus increment to keep track of insertions
        if verb:
            print 'sizel: inserting ' + str(dupii) + ' before index ' + str(dupii+1+incr-unqidx[0])
        sizel.insert(dupii+1+incr-unqidx[0],dupii)

        # dupesinsert are indices into raw array
        if verb:
            print 'heightl: inserting ' + str(heightraw[dupesinsert[dii]+1]) + \
                ' before index ' + str(dupii+1+incr-unqidx[0]) +\
                ' (dupesinsert= ' + str(dupesinsert[dii]+1) + ')'
        heightl.insert(dupii+1+incr-unqidx[0], heightraw[dupesinsert[dii]+1])
        sizerawl.insert(dupii+1+incr-unqidx[0], sizeraw[dupesinsert[dii]+1]) 
        notesl[dupii+incr-unqidx[0]] = 'DUPLICATE SIZE'
        notesl.insert(dupii+1+incr-unqidx[0], 'DUPLICATE SIZE')

        incr+=1 # now we have one additional element in array. increase insertion index
        dii+=1 # index into dupesinsert

    # convert back to numpy arrays
    sizeary = np.array(sizel)
    heightary = np.array(heightl)
    sizerawary = np.array(sizerawl)                   
    notesary = np.array(notesl)

    if plotfig:

        fig = plt.figure(figsize=(9,5)) 
        plt.plot(sizeraw,seldf[heightcol],marker='s',color='r',linestyle='none')
        plt.plot(sizearyunq,heightaryunq,marker='o',markersize=6,fillstyle='none')        
        plt.plot(sizeary,heightary,marker='^',color='g',linestyle='none')
        plt.xlim((sizelims[0]-2,sizelims[1]+2))
        plt.title(fname)
        plt.xlabel(selcol)
        plt.ylabel(heightcol)
        for dii,dup in enumerate(dupes):
            plt.axvline(x=dup,color='k',linestyle='--')
            if verb:
                print 'duplicate ' + selcol + ', ' + heightcol + ': (' + str(dup) + ',' +\
                    str(heightary[dup+dii-sizelims[0]]) +\
                    '), (' + str(dup) + ',' + str(heightary[dup+1+dii-sizelims[0]]) + ')'

        plt.legend(('Raw','Processed: No duplicates',
                    'Processed: w/ Duplicates','Duplicates'), fancybox=True,
                   framealpha=0.5, frameon=False)
        if printtofile:
            splits=fname.split('/')[-1].split('.')
            if verb:
                print 'Printing figure: peakal_' + splits[0] + '.pdf'
            fig.savefig('peakal_' + splits[0] + '.pdf')

    return ((sizeary,heightary),(sizearyunq,heightaryunq),(sizeraw,heightraw), sizerawary, notesary)


def plot_sample_diff(samp1, samp2, samp1unq, samp2unq, sampnames=None):

    if sampnames == None:
        sampnames=('1','2')

    fig,axs=plt.subplots(2,1)
    fig.set_size_inches((9,5))
    ax=axs[0]
    ax.plot(samp1[0],samp1[1],color='r',marker='s')
    ax.plot(samp2[0],samp2[1],color='b',marker='.')
    ax.set_ylabel('Height')
    ax.legend((sampnames[0],sampnames[1]),frameon=False,loc='best')
    ax.grid(True)
    ax=axs[1]
    ax.plot(samp1unq[0],samp1unq[1]-samp2unq[1],color='k',marker='.')
    ax.set_ylabel('Height')
    ax.legend((sampnames[0] + ' - ' + sampnames[1],),frameon=False,loc='best')
    ax.axhline(y=0,color='k',linestyle='--')
    ax.grid(True)

    return fig,axs


def read_samples(fname, sampcol='Sample File Name', delimiter='\t'):
    """ returns a list of sample dataframes
           
    """
    datdf=pd.read_csv(fname,delimiter=delimiter)
    cols = datdf.keys()
    sampkeys = np.unique(datdf[sampcol].values)
    
    samples = []
    for sampkey in sampkeys:
        sampdf=pd.DataFrame(datdf.values[np.array(datdf[sampcol]==sampkey)],
                       columns=cols)
        samples.append(sampdf)

    return samples

def get_sample_name(sample, selcol='Sample File Name'):
    
    return sample[selcol].values[0]



def test_main(fname, sizerange, sampcol='Sample File Name', sizecol='Size', 
              heightcol='Height', delimiter='\t', plotfig=False, printtofile=False, 
              verb=False):
    """
        read samples from fname, within the given sizerange (integers, inclusive of endpoints)

            Expecting 1-2 samples within file
    """

    if verb:
        print 'Reading samples'
    samples = read_samples(fname, sampcol=sampcol, delimiter=delimiter)
    
    if type(sizerange) == str:
        splits=sizerange.split(',')
        sizelims=np.array([np.int(splits[0]), np.int(splits[1])])
    else:
        sizelims=sizerange

    # multiple samples within in each file:
    # proc is [#samples][size/height]
    proc={}; procun={}; procraw={}; sampnames={}; sizerawdt={}; notesdt={}
    colstack = []
    for sii,samp in enumerate(samples):
        sampnames[sii]=get_sample_name(samp)
        if verb:
            print '== PROCESSING SAMPLE: ' + sampnames[sii]
        proc[sii],procun[sii],procraw[sii], sizerawdt[sii], notesdt[sii] = process_sample(samp, sizelims, selcol=sizecol,
                                                                                          heightcol=heightcol,
                                                                                          delimiter=delimiter,
                                                                                          plotfig=plotfig,
                                                                                          verb=verb,
                                                                                          printtofile=printtofile)

        # stack up the output columns
        colstack.append(np.vstack(((sampnames[sii],)*len(proc[sii][0]),
                                   proc[sii][0], sizerawdt[sii], 
                                   proc[sii][1], notesdt[sii])))

    outstack = np.hstack(colstack)
    outdf = pd.DataFrame(outstack.T,columns=(sampcol,sizecol,sizecol + ' (raw)', heightcol, 'Notes'))
    outdf.to_csv(fname+'_out'+ext,sep=delimiter)

    if plotfig and len(samples)==2:
        # @@@ plots just the first 2 samples
        fig,axs = plot_sample_diff(proc[0], proc[1], procun[0], procun[1], sampnames=sampnames)



printtofile=False

basepath='./Ryan/'
ext = '.txt'
fnameb = basepath + 'B_peak_height'
fnamea = basepath + 'A_peak_height'
fname1 = basepath + 'set1'
fname2 = basepath + 'set2'
fname3 = basepath + 'set3'
fnameh = basepath + 'hnsdata'
fnamebug = basepath + 'hns_c_03_test'
delimiter='\t'

fin = fnamebug+ext
#sizelims=[100,200]
sizelims='100,200'
test_main(fin, sizelims, plotfig=True, verb=True)

# in HNS file:
# check C_C03.fsa
# sizes ~ 150-156
# for some reason, duplicates found is incorrect:
#   Duplicates found: [122 156 200 150 197 200]

"""
samples1 = read_samples(fname1+ext)

# multiple samples within in each file:
# proc is [#samples][size/height]
# FILE1
proc={}; procun={}; procraw={}; sampnames={}; sizerawdt={}; notesdt={}
for sii,samp in enumerate(samples1):
    sampnames[sii]=get_sample_name(samp)
    proc[sii],procun[sii],procraw[sii], sizerawdt[sii], notesdt[sii] = process_sample(samp,sizelims,printtofile=printtofile)
fig,axs = plot_sample_diff(proc[0], proc[1], procun[0], procun[1], sampnames=sampnames)

out1s1 = np.vstack(((sampnames[0],)*len(proc[0][0]),proc[0][0],sizerawdt[0],proc[0][1],notesdt[0]))
out1s2 = np.vstack(((sampnames[1],)*len(proc[1][0]),proc[1][0],sizerawdt[1],proc[1][1],notesdt[1]))
out1 = np.hstack((out1s1,out1s2))
out1df = pd.DataFrame(out1.T,columns=('Sample File Name','Size','Size (raw)', 'Height', 'Notes'))
out1df.to_csv(fname1+'_out'+ext,sep=delimiter)


# FILE2
samples2 = read_samples(fname2+ext)
proc2={}; procun2={}; procraw2={}; sampnames2={}; sizerawdt2={}; notesdt2={}
for sii,samp in enumerate(samples2):
    sampnames2[sii]=get_sample_name(samp)
    proc2[sii],procun2[sii],procraw2[sii], sizerawdt2[sii], notesdt2[sii]  = process_sample(samp,sizelims,printtofile=printtofile)
plot_sample_diff(proc2[0], proc2[1], procun2[0], procun2[1], sampnames=sampnames2)    

# FILE3
samples3 = read_samples(fname3+ext)
proc3={}; procun3={}; procraw3={}; sampnames3={}; sizerawdt3={}; notesdt3={}
for sii,samp in enumerate(samples3):
    sampnames3[sii]=get_sample_name(samp)
    proc3[sii],procun3[sii],procraw3[sii], sizerawdt3[sii], notesdt3[sii] = process_sample(samp,sizelims,printtofile=printtofile)
plot_sample_diff(proc3[0], proc3[1], procun3[0], procun3[1], sampnames=sampnames3)    


samplesa = read_samples(fnamea+ext)
samplesb = read_samples(fnameb+ext)

(proca, procuna, procrawa, sizerawdta, notesdta) = process_sample(samplesa[0],sizelims,printtofile=printtofile)
(procb, procunb, procrawb, sizerawdtb, notesdtb) = process_sample(samplesb[0],sizelims,printtofile=printtofile)





fig,axs=plt.subplots(2,1)
fig.set_size_inches((9,5))
ax=axs[0]
ax.plot(proca[0],proca[1],color='r',marker='s')
ax.plot(procb[0],procb[1],color='b',marker='.')
ax.set_ylabel('Height')
ax.legend(('A','B'),frameon=False,loc='best')
ax.grid(True)
ax=axs[1]
ax.plot(procuna[0],procuna[1]-procunb[1],color='k',marker='.')
ax.set_ylabel('Height')
ax.legend(('A-B',),frameon=False,loc='best')
ax.axhline(y=0,color='k',linestyle='--')
ax.grid(True)
#ax.set_title('A-B')
if printtofile:
    fig.savefig('peakal_AB_A-B.pdf')

"""

"""
sizerange=np.arange(sizelims[0],sizelims[1]+1) # @@@@ NEW question: range inclusive of endpoints?

datadf = pd.read_csv(fnamea,delimiter="\t")
datbdf = pd.read_csv(fnameb,delimiter="\t")


#Aret = datdf.values[datdf['Sample File Name'].values == 'A_A03.fsa']
#Bret = datdf.values[datdf['Sample File Name'].values == 'B_B03.fsa']

# A DATASET: ==================
acols=datadf.keys()

#adat = pd.DataFrame(datadf.values[np.logical_and(datadf['Size'] >= sizelims[0], 
#                    datadf['Size'] <= sizelims[1])], columns=acols)
adat = datadf
asizerawall = datadf['Size'] # save original size
# save raw selection
asizeraw = asizerawall[np.logical_and(asizerawall.round()>=sizelims[0],
                                      asizerawall.round()<=sizelims[1])] 
# Round size data
adat['Size']=datadf['Size'].round()
# Select data based on user Size range:
selbool=np.array(np.logical_and(adat['Size'] >= sizelims[0],adat['Size'] <= sizelims[1]))
asel = pd.DataFrame(adat.values[selbool], columns=acols)

# Convert Size to integer vals and use as an index
arawidx=asel['Size'].values.astype(int) 
anumsize=len(arawidx)

# Find duplicates:
# Use dupe indices to insert data later (this only works b/c arawidx is sorted)
dupes=arawidx[arawidx[1:] == arawidx[:-1]]

# remove duplicate sizes (indices)
# # reconidx reconstructs original, retidx are the indices that result in unique array
aidxunq, aretidx, areconidx = np.unique(arawidx,return_index=True,return_inverse=True) 
# use retidx to set height data associatd w/ unique size indices

dupesinsert=np.zeros((len(dupes)))
for dii,dup in enumerate(dupes):
    dupesinsert[dii] = np.int(cutl.find_nearest(arawidx,dup))
#dupesinsert=areconidx[areconidx[1:] == areconidx[:-1]] #@@ no this doesn't work for my purposes!


# Get Heights that corresponds to Sizes
aheightraw=asel['Height'].values


# Initialize final height array and final size array 
aheightary=np.zeros(len(sizerange)) 
asizeary=np.arange(sizelims[0],sizelims[1]+1)

aheightraw = aheightraw.astype(type(aheightraw[0]))
aheightary[aidxunq-sizelims[0]] = aheightraw[aretidx]


plt.figure() 
plt.plot(asizeary,aheightary,marker='o')
plt.plot(asizeraw.values,asel['Height'],marker='s',color='r')
plt.legend(('Processed: No duplicates','Raw'))

# ADD DUPLICATES BACK IN:
aheightl = list(aheightary)
asizel = list(asizeary)
dii=0
for dupii in dupes: # dupe size values (also indices into full size array)

    # insert the value before the index
    # list.insert(index,value)
    print 'sizel: inserting ' + str(dupii) + ' before index ' + str(dupii+2-aidxunq[0])
    asizel.insert(dupii+2-aidxunq[0],dupii)

    # dupesinsert are indices into raw array
    print 'heightl: inserting ' + str(aheightraw[dupesinsert[dii]+1]) + \
        ' before index ' + str(dupii+2-aidxunq[0]) + ' (dupesinsert= ' + str(dupesinsert[dii]+1) + ')'
    aheightl.insert(dupii+2-aidxunq[0], aheightraw[dupesinsert[dii]+1])

    dii+=1 # index into dupesinsert

asizearydup=np.array(asizel)
aheightarydup=np.array(aheightl)

plt.figure() 
plt.plot(asizeary,aheightary,marker='o',fillstyle='none')
plt.plot(asizeraw.values,asel['Height'],marker='s',color='r')
plt.plot(asizearydup,aheightarydup,marker='*',color='g',linestyle='none')
plt.xlim((sizelims[0]-2,sizelims[1]+2))
for dup in dupes:
    plt.axvline(x=dup,color='k',linestyle='--')
plt.legend(('Processed: No duplicates','Raw','Processed: with dupes','Dupes'))








# B DATASET: =================
bcols=datbdf.keys()
#bdat = pd.DataFrame(datbdf.values[np.logical_and(datbdf['Size'] >= sizerange[0], datbdf['Size'] <= sizerange[1])], columns=bcols)

bdat = datbdf
bdat['Size']=datbdf['Size'].round()
bsel = pd.DataFrame(bdat.values[np.logical_and(bdat['Size'] >= sizerange[0], bdat['Size'] <= sizerange[1])], columns=bcols)

"""


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


"""
def myround(vals):

    rndvals=np.zeros(len(vals))
    for ii in np.arange(0,len(vals)):
        rndvals[ii] = np.round(vals[ii])
        #print vals[ii],rndvals[ii]

    return rndvals



def fillarray(size, height,sizerange, verb=True):

    # take size and height arrays, round Size. 
    #    If an index in sizerange is missing in size, enter 0 for height
    

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

"""
