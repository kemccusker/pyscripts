"""
    hiatus.py
         2/27/2014: plot of dummy trend values as example
         
"""

import matplotlib.pyplot as plt # for basic plotting
import matplotlib.cm as cm

plt.close("all")
plt.ion()
#         ANN, cold, warm  (C/dec)
glob15 = [0.05, 0.00, 0.09]  # GST
glob15ci = [0.08, 0.10, 0.08] # symmetric confidence intervals
glob15sig = [0, 0, 1]

land15 = [0.09, 0.02, 0.18] # GSAT
land15ci = [0.13,0.19,0.14]
land15sig = [0, 0, 1]

ocean15 = [0.03,-0.02,0.05]
ocean15ci = [0.07,0.08,0.07]
ocean15sig = [0, 0, 0]

fig, ax1 = plt.subplots(figsize=(10,6))
plt.plot([0,4],[0,0],color='k')

plt.plot([1,1],[glob15[0]-glob15ci[0],glob15[0]+glob15ci[0]],color='gray')
plt.plot([1.1,1.1],[glob15[1]-glob15ci[1],glob15[1]+glob15ci[1]],color='gray')
plt.plot([1.2,1.2],[glob15[2]-glob15ci[2],glob15[2]+glob15ci[2]],color='gray')
plt.plot([1,1.1],glob15[0:2],linestyle='none',marker='s',color='k',markersize=10,mew=2,mfc='white',fillstyle='none')#,fillstyle=('none','none','full'))
plt.plot([1.2],glob15[2],linestyle='none',marker='s',color='k',markersize=10,mew=2)

plt.plot([2,2],[land15[0]-land15ci[0],land15[0]+land15ci[0]],color='gray')
plt.plot([2.1,2.1],[land15[1]-land15ci[1],land15[1]+land15ci[1]],color='gray')
plt.plot([2.2,2.2],[land15[2]-land15ci[2],land15[2]+land15ci[2]],color='gray')
plt.plot([2,2.1],land15[0:2],linestyle='none',marker='s',color='brown',markersize=10,mew=2,mfc='white',fillstyle='none')
plt.plot([2.2],land15[2],linestyle='none',marker='s',color='brown',mec='brown',markersize=10,mew=2) #markeredgecolor

plt.plot([3,3],[ocean15[0]-ocean15ci[0],ocean15[0]+ocean15ci[0]],color='gray')
plt.plot([3.1,3.1],[ocean15[1]-ocean15ci[1],ocean15[1]+ocean15ci[1]],color='gray')
plt.plot([3.2,3.2],[ocean15[2]-ocean15ci[2],ocean15[2]+ocean15ci[2]],color='gray')
plt.plot([3,3.1,3.2],ocean15,linestyle='none',marker='s',color='blue',markersize=10,mew=2,mfc='white',fillstyle='none')#markeredgewidth=mew

plt.title('1998-2012 Trends (C). Closed marker=significant at 90%. Uncertainty interval: 5-95%')
plt.ylim([-.2,0.4])
plt.xlim([0,4])
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)
#ax1.set_xticks([1,2,3])
#ax1.set_xticklabels(('GST','GSAT','GSST'))# put these on top of plot
ax1.set_xticks([1,1.1,1.2,2,2.1,2.2,3,3.1,3.2])
ax1.set_xticklabels(('ann','cold','warm','ann','cold','warm','ann','cold','warm'),rotation=45)

# angle xticklabels for seasons ANN, cold warm
## ax1.text(1,0.37,'ann',rotation=45)
## ax1.text(1.1,0.37,'cold',rotation=45)
## ax1.text(1.2,0.37,'warm',rotation=45)
ax1.text(1,0.37,'global') #GST
ax1.text(2,0.37,'land') # GSAT
ax1.text(3,0.37,'ocean') #GSST

plt.savefig('hiatustrends.pdf')


# partial example from: http://matplotlib.org/examples/pylab_examples/boxplot_demo2.html
#
## fig, ax1 = plt.subplots(figsize=(10,6))
## fig.canvas.set_window_title('A Boxplot Example')
## plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

## bp = plt.boxplot(data, notch=0, sym='+', vert=1, whis=1.5)
## plt.setp(bp['boxes'], color='black')
## plt.setp(bp['whiskers'], color='black')
## plt.setp(bp['fliers'], color='red', marker='+')

## # Add a horizontal grid to the plot, but make it very light in color
## # so we can use it for reading data values but not be distracting
## ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
##               alpha=0.5)

## # Hide these grid behind plot objects
## ax1.set_axisbelow(True)
## ax1.set_title('Comparison of IID Bootstrap Resampling Across Five Distributions')
## ax1.set_xlabel('Distribution')
## ax1.set_ylabel('Value')
