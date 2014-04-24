"""
cccmacmaps.py
  My own discrete colormaps

2/27/2014

"""
import numpy as np # for array handling
import matplotlib.colors as col
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt # for basic plotting


def get_cccmacm():

    # do an if elif elif else for colormaps
    pass



def list_cccmacms():
    # print the colormap name options
    return ['kem_w20','blue2red_w20','blue2red_20', 'red2blue_w20', 'turq2orange_16',\
            'brown2blue_12', 'brown2blue_12g', 'brown2blue_12w',\
            'brown2blue_16w', 'brown2blue_16g', \
            'blue2blue_9', 'blue2blue_w10', 'blue2blue_bw10']


def show_cccmacms():
    # plot the colormaps as the demo example
    # http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps

    register_cccmacms()
    names = list_cccmacms()
    
    # COPIED BELOW CODE
    # base code from http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps
    matplotlib.rc('text', usetex=False)
    a=np.outer(np.arange(0,1,0.01),np.ones(10))   # pseudo image data
    f=plt.figure(figsize=(7,9))
    f.subplots_adjust(top=0.8,bottom=0.05,left=0.01,right=0.99)
    # get list of all colormap names
    # this only obtains names of built-in colormaps:
    maps=[m for m in cm.datad if not m.endswith("_r")]
    # use undocumented cmap_d dictionary instead
    maps = [m for m in cm.cmap_d if not m.endswith("_r")]
    maps.sort()
    # determine number of subplots to make
    l=len(maps)+1
    if names is not None: l=len(names)  # assume all names are correct!
    # loop over maps and plot the selected ones
    i=0
    for m in maps:
        if names is None or m in names:
            i+=1
            ax = plt.subplot(1,l,i)
            ax.axis("off")
            plt.imshow(a,aspect='auto',cmap=cm.get_cmap(m),origin="lower")
            plt.title(m,rotation=90,fontsize=10,verticalalignment='bottom')
            
    # END COPIED CODE


def register_cccmacms(cmap='all'):
    """create my personal colormaps with discrete colors and register them.
       default is to register all of them. can also specify which one.
       (@@ input arg cmap not implemented yet 2/27/14)
    """
    #print 'registering cmaps'
    
    # define individual colors as RGB triples
    # from colorwheel.m
    # =============================================
    # kem_w20 (20)  OR blue2red_w20
    # blueish at top, white in middle, reddish at bottom

    cpool = np.array([ [10,50,120], \
                       [15,75,165], \
                       [30,110,200],\
                       [60,160,240],\
                       [80,180,250],\
                       [130, 210, 255],\
                       [160, 230, 255],\
                       [190, 235, 255],\
                       [210, 245, 255],\
                       [255, 255, 255],\
                       
                       [255, 255, 255],\
                       [250, 240, 150],\
                       [255, 222, 100],\
                       [255, 192, 60], \
                       [255, 160, 0], \
                       [255, 96, 0], \
                       [255, 50, 0], \
                       [225, 20, 0], \
                       [192, 0, 0], \
                       [165, 0, 0]],\
                       dtype=float)
    
    kem_w20 = (cpool/255.)
    thecmap = col.ListedColormap(kem_w20,'kem_w20')
    cm.register_cmap(cmap=thecmap)

    # really should call this blue2red_w20
    blue2red_w20 = (cpool/255.)
    thecmap = col.ListedColormap(blue2red_w20,'blue2red_w20')
    cm.register_cmap(cmap=thecmap)

    # ============================================
    # red2blue_w20 (the above, flipped)
    #
    red2blue_w20 = np.flipud(cpool)/255.
    thecmap = col.ListedColormap(red2blue_w20,'red2blue_w20')
    cm.register_cmap(cmap=thecmap)
    

    # =============================================
    # blue2red_20
    # blueish at top, reddish at bottom

    cpool = np.array([ [10,0,110], \
                       [20,50,130], \
                       [15,75,165], \
                       [30,110,200],\
                       [60,160,240],\
                       [80,180,250],\
                       [130, 210, 255],\
                       [160, 230, 255],\
                       [190, 235, 255],\
                       #[210, 245, 255],\
                       [215, 250, 255],\
                       
                       [255, 250, 220],\
                      # [250, 240, 150],\
                       [255, 222, 100],\
                       [255, 192, 60], \
                       [255, 160, 0], \
                       [255, 96, 0], \
                       [255, 50, 0], \
                       [225, 20, 0], \
                       [192, 0, 0], \
                       [165, 10, 0],\
                       [110, 0, 0]],\
                       dtype=float)
    
    blue2red_20 = (cpool/255)
    thecmap = col.ListedColormap(blue2red_20,'blue2red_20')
    cm.register_cmap(cmap=thecmap)


    # =============================================
    # turq2orange_16 (formerly named blue2orange_18)
    #   seagreen/turquoise at top, rusty orange at bottom
    cpool = np.array([ [0, 102, 102], \
                       [0, 153, 153], \
                       [0, 204, 204], \
#                       [0, 255, 255], \
#                       [51, 255, 255], \
                       [101, 255, 255], \
                       [153, 255, 255], \
                       [178, 255, 255], \
                       [203, 255, 255], \
                       
                       [229, 255, 255], \
                       [255, 229, 203], \
                       [255, 202, 153], \
                       [255, 173, 101], \
                       [255, 142, 51], \
                       [255, 110, 0], \
                       [204, 85, 0], \
                       [153, 61, 0], \
                       [102, 39, 0] ],\
                     dtype=float)

    turq2orange_16 = cpool/255
    thecmap = col.ListedColormap(turq2orange_16,'turq2orange_16')
    cm.register_cmap(cmap=thecmap)

    # =============================================
    # brown2blue_12 (formerly named kembr2bl)
    #    brown at top, blue at bottom
    cpool = np.array([ [51, 25, 0], \
                       [102, 47, 0],\
                       [153, 96, 53], \
                       [204, 155, 122], \
                       [216, 175, 151], \
                       [242, 218, 205], \

                       [224, 255, 255], \
                       [170, 247, 255], \
                       [114, 217, 255], \
                       [63, 160, 255], \
                       [38, 77, 255], \
                       [41, 10, 216] ], \
                     dtype=float)

    brown2blue_12 = cpool/255
    thecmap = col.ListedColormap(brown2blue_12,'brown2blue_12')
    cm.register_cmap(cmap=thecmap)

    # =============================================
    # brown2blue_12g (formerly named kemgray)
    #    Same as brown2blue_12 but gray in middle
    #    brown at top, gray in middle, blue at bottom
    cpool = np.array([ [51, 25, 0], \
                       [102, 47, 0],\
                       #[153, 96, 53], \
                       [204, 155, 122], \
                       [216, 175, 151], \
                       [242, 218, 205], \
                       
                       [240, 240, 240],\
                       [240, 240, 240],\

                       [224, 255, 255], \
                       [170, 247, 255], \
                       [114, 217, 255], \
                      # [63, 160, 255], \
                       [38, 77, 255], \
                       [41, 10, 216] ], \
                     dtype=float)

    brown2blue_12g = cpool/255
    thecmap = col.ListedColormap(brown2blue_12g,'brown2blue_12g')
    cm.register_cmap(cmap=thecmap)

    # =============================================
    # brown2blue_12w
    #    Same as brown2blue_12g but white in middle
    #    brown at top, white in middle, blue at bottom
    cpool = np.array([ [51, 25, 0], \
                       [102, 47, 0],\
                       #[153, 96, 53], \
                       [204, 155, 122], \
                       [216, 175, 151], \
                       [242, 218, 205], \

                       [255, 255, 255],\
                       [255, 255, 255],\

                       [224, 255, 255], \
                       [170, 247, 255], \
                       [114, 217, 255], \
                      # [63, 160, 255], \
                       [38, 77, 255], \
                       [41, 10, 216] ], \
                     dtype=float)

    brown2blue_12w = cpool/255
    thecmap = col.ListedColormap(brown2blue_12w,'brown2blue_12w')
    cm.register_cmap(cmap=thecmap)

# =============================================
    # brown2blue_16w
    #    New browns. Similar blues to brown2blue_12w. white in middle
    #    brown at top, white in middle, blue at bottom
    cpool = np.array([ [51, 25, 0],\
                       #[72, 60, 50], \
                       #[128, 70, 27], \
                       [100, 68, 35],\
                       [130, 102, 68], \
                       [159, 139, 112], \
                       [195, 176, 145], \
                      # [193, 154, 107], \
                      # [230, 210, 200],\
                       [215, 200, 181], \
                       [232, 228, 215], \

                       [255, 255, 255],\
                       [255, 255, 255],\

                       [224, 255, 255], \
                       [170, 247, 255], \
                       [114, 217, 255], \
                       [63, 160, 255], \
                       [38, 77, 255], \
                       [41, 10, 216], \
                       #[20, 50, 130], \
                       [10, 0, 110] ],\
                     dtype=float)
# [10,0,110], \ # darkest blues from blue2red_20
# [20,50,130]
                     
# http://en.wikipedia.org/wiki/Shades_of_brown
# 72, 60, 50 # taupe. very dark.
# 107, 68, 35 # brown-nose
# 130, 102, 68 # raw umber
# 128, 70, 27 # russett
# 159, 139, 112 # beaver
# 195, 176, 145 # khaki
#  [255, 222, 173], # navajo white http://www.rapidtables.com/web/color/RGB_Color.htm
# 193, 154, 107 wood brown

    brown2blue_16w = cpool/255
    thecmap = col.ListedColormap(brown2blue_16w,'brown2blue_16w')
    cm.register_cmap(cmap=thecmap)   

# =============================================
    # brown2blue_16g
    #    Same as brown2blue_16w but gray middle
    #    brown at top, gray in middle, blue at bottom
    cpool = np.array([ [51, 25, 0],\
                       #[72, 60, 50], \
                       #[128, 70, 27], \
                       [100, 68, 35],\
                       [130, 102, 68], \
                       [159, 139, 112], \
                       [195, 176, 145], \
                      # [193, 154, 107], \
                      # [230, 210, 200],\
                       [215, 200, 181], \
                       [232, 228, 215], \

                       [240, 240, 240],\
                       [240, 240, 240],\

                       [224, 255, 255], \
                       [170, 247, 255], \
                       [114, 217, 255], \
                       [63, 160, 255], \
                       [38, 77, 255], \
                       [41, 10, 216], \
                       #[20, 50, 130], \
                       [10, 0, 110] ],\
                     dtype=float)
    brown2blue_16g = cpool/255
    thecmap = col.ListedColormap(brown2blue_16g,'brown2blue_16g')
    cm.register_cmap(cmap=thecmap) 


    # ===============================================
    # blue2blue_9
    #     light to dark blue
    cpool = np.array([ [229, 255, 255], \
                       [204, 250, 255], \
                       [178, 242, 255], \
                       [153, 229, 255], \
                       [127, 212, 255], \
                       [101, 191, 255], \
                       [76, 165, 255], \
                       [50, 101, 255], \
                       [0, 63, 255] ], \
                     dtype=float)

    blue2blue_9 = cpool/255
    thecmap = col.ListedColormap(blue2blue_9,'blue2blue_9')
    cm.register_cmap(cmap=thecmap)

    # ===============================================
    # blue2blue_w10
    #     Same as blue2blue_9, but first color is white
    #     white, light to dark blue
    cpool = np.array([ [255, 255, 255], \
                       [229, 255, 255], \
                       [204, 250, 255], \
                       [178, 242, 255], \
                       [153, 229, 255], \
                       [127, 212, 255], \
                       [101, 191, 255], \
                       [76, 165, 255], \
                       [50, 101, 255], \
                       [0, 63, 255] ], \
                     dtype=float)

    blue2blue_w10 = cpool/255
    thecmap = col.ListedColormap(blue2blue_w10,'blue2blue_w10')
    cm.register_cmap(cmap=thecmap)    

    # ===============================================
    # blue2blue_bw10
    #     Similar to blue2blue_w10, but 1st color is light brown, then white
    #     light brown, white, light to dark blue
    cpool = np.array([ [242, 218, 205], \
                       [255, 255, 255], \
                       [200, 255, 255], \
#                       [229, 255, 255], \
#                       [204, 250, 255], \
                       [178, 242, 255], \
                       [153, 229, 255], \
                       [127, 212, 255], \
                       [101, 191, 255], \
                       [76, 165, 255], \
                       [50, 101, 255], \
                       [0, 63, 255] ], \
                     dtype=float)

    blue2blue_bw10 = cpool/255
    thecmap = col.ListedColormap(blue2blue_bw10,'blue2blue_bw10')
    cm.register_cmap(cmap=thecmap)  
                      


# this means, if we are running this module as main script
#    (not importing from another), register and show
if __name__ == "__main__":
    register_cccmacms()
    show_cccmacms()
else:
    register_cccmacms()




## # # define the colormap function # # #
## def discrete_cmap(N=22):
##     """create a colormap with N (N<15) discrete colors and register it"""

##     # define individual colors as RGB triples
##     # from colorwheel.m
##     # kem_w22 (22)
##     # blueish at top, white in middle, reddish at bottom

##     cpool = np.array([ [10,50,120], \
##                        [15,75,165], \
##                        [30,110,200],\
##                        [60,160,240],\
##                        [80,180,250],\
##                        [130, 210, 255],\
##                        [160, 230, 255],\
##                        [190, 235, 255],\
##                        [210, 245, 255],\
##                        [255, 255, 255],\
##                        [255, 255, 255],\
##                        [250, 240, 150],\
##                        [255, 222, 100],\
##                        [255, 192, 60], \
##                        [255, 160, 0], \
##                        [255, 96, 0], \
##                        [255, 50, 0], \
##                        [225, 20, 0], \
##                        [192, 0, 0], \
##                        [165, 0, 0]],\
##                        dtype=float)
    
##     kem_w22 = (cpool/255);
##     thecmap = col.ListedColormap(kem_w22,'kem_w22')
##     cm.register_cmap(cmap=thecmap)
