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
    return ['kem_w20','blue2red_w20','blue2red_20', \
            'blue2red_w10','blue2red_w11','red2blue_w20',\
            'red2blue_w20sic', 'red2blue_20',\
            'turq2orange_16',\
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
    

    # ============================================
    # ======= specialty red2blue for SIC anomalies
    # red2blue_w20sic
    #       lots of white in middle
    
    # start with blue2red as above
    # anomalies should go from -30, 30 (if %) for white from -10 to 10
    # anoms from -25 to 25 for white from +-5
    cpool = np.array([ [10,50,120], \
                       [15,75,165], \
                       [30,110,200],\
                       [60,160,240],\
                       [80,180,250],\
                       #[130, 210, 255],\
                       [160, 230, 255],\
                       #[190, 235, 255],\

                       [210, 245, 255],\
                       [210, 245, 255],\

                       #[255, 255, 255],\
                       #[255, 255, 255],\
                       [255, 255, 255],\
                       [255, 255, 255],\
                       
                       [255, 255, 255],\
                       [255, 255, 255],\
                       #[255, 255, 255],\
                       #[255, 255, 255],\

                       [250, 240, 150],\
                       [250, 240, 150],\

                       [255, 222, 100],\
                       #[255, 192, 60], \
                       [255, 160, 0], \
                       [255, 96, 0], \
                       [255, 50, 0], \
                       [225, 20, 0], \
                       #[192, 0, 0], \
                       [165, 0, 0]],\
                       dtype=float)
    red2blue_w20sic = np.flipud(cpool)/255.
    thecmap = col.ListedColormap(red2blue_w20sic,'red2blue_w20sic')
    cm.register_cmap(cmap=thecmap)


    # =============================================
    # blue2red_w10
    # blueish at top, white in middle, reddish at bottom

    cpool = np.array([ #[10,50,120], \
                       [15,75,165], \
                       #[30,110,200],\
                       [60,160,240],\
                       #[80,180,250],\
                       [130, 210, 255],\
                       #[160, 230, 255],\
                       [190, 235, 255],\
                       #[210, 245, 255],\
                       [255, 255, 255],\
                       
                       [255, 255, 255],\
                       #[250, 240, 150],\
                       [255, 222, 100],\
                       #[255, 192, 60], \
                       [255, 160, 0], \
                       #[255, 96, 0], \
                       [255, 50, 0], \
                       #[225, 20, 0], \
                       [192, 0, 0]], \
                       #[165, 0, 0]],\
                       dtype=float)
    
    #  blue2red_w10
    blue2red_w10 = (cpool/255.)
    thecmap = col.ListedColormap(blue2red_w10,'blue2red_w10')
    cm.register_cmap(cmap=thecmap)

    # =============================================
    # blue2red_w11
    # blueish at top, white in middle, reddish at bottom

    cpool = np.array([ [10,50,120], \
                       #[15,75,165], \
                       [30,110,200],\
                       #[60,160,240],\
                       [80,180,250],\
                       #[130, 210, 255],\
                       [160, 230, 255],\
                       #[190, 235, 255],\
                       [210, 245, 255],\

                       [255, 255, 255],\

                       [250, 240, 150],\
                       #[255, 222, 100],\
                       [255, 192, 60], \
                       #[255, 160, 0], \
                       [255, 96, 0], \
                       #[255, 50, 0], \
                       [225, 20, 0], \
                       #[192, 0, 0], \
                       [165, 0, 0]],\
                       dtype=float)
    
    #  blue2red_w11
    blue2red_w11 = (cpool/255.)
    thecmap = col.ListedColormap(blue2red_w11,'blue2red_w11')
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

    # ============================================
    # red2blue_20 (the above, flipped)
    #
    red2blue_20 = np.flipud(cpool)/255.
    thecmap = col.ListedColormap(red2blue_20,'red2blue_20')
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
                      

def get_linecolor(lckey):
    """ get_linecolor(lckey)
             return the RGB array defining a color

                  'firebrick1': np.array([255, 48, 48])/255.,
                  'firebrick': np.array([178, 34, 34])/255.,
                  'orangered': np.array([255, 69, 0])/255.,
                  'orangered4': np.array([139, 37, 0])/255.,
                  'darkgoldenrod1': np.array([255, 185, 15])/255.,
                  'darkgoldenrod': np.array([184, 134, 11])/255.,
                  'darkolivegreen1': np.array([202, 255, 112])/255.,
                  'darkolivegreen3': np.array([162, 205, 90])/255.,
                  'darkseagreen': np.array([143, 188, 143])/255.,
                  'darkseagreen4': np.array([105, 139, 105])/255.,
                  'mediumpurple1': np.array([171, 130, 255])/255.,
                  'mediumpurple4': np.array([93, 71, 139])/255.,
                  'dodgerblue': np.array([30, 144, 255])/255.,
                  'mediumblue': np.array([0, 0, 205])/255.,
                  'limegreen': np.array([50, 205, 50])/255., 
                  'lightsteelblue3': np.array([162, 181, 205])/255., # more grey looking
                  'lightsteelblue4': np.array([110, 123, 139])/255., # more grey looking
                  'steelblue3': np.array([79, 148, 205])/255.,  # more blue looking
                  'steelblue4': np.array([54, 100, 139])/255.,
                  'pink2': np.array([238, 169, 184])/255.,
                  'pink4': np.array([139, 99, 108])/255.,
                  'chocolate1': np.array([255, 127, 36])/255.,
                  'chocolate4': np.array([139, 69, 19])/255.
                  'barnred': np.array([124, 10, 2]/255.,
                  'yelloworange': np.array([255, 155, 0])/255.,
                  'darkyellow': np.array([242, 228, 40])/255.,
                  'magenta': np.array([112, 0, 102])/255.,
                  'violet': np.array([85, 0, 140])/255.,
                  'deepskyblue': np.array([0,191,255])/255.,
                  'skyblue': np.array([135, 206, 250])/255.,
                  'turq': np.array([0, 153, 153])/255.
                  'niceblue': np.array([0,102,204])/255.

                   http://web.njit.edu/~kevin/rgb.txt.html
                   http://en.wikipedia.org/wiki/Shades_of_red
                   http://www.rapidtables.com/web/color/RGB_Color.htm

    """
    linecolors = get_linecolorwheel()
    
    return linecolors[lckey]


def get_linecolorwheel():
    """ get_linecolorwheel()

            return the dictionary of all linecolors

             http://web.njit.edu/~kevin/rgb.txt.html
             http://en.wikipedia.org/wiki/Shades_of_red
             http://www.rapidtables.com/web/color/RGB_Color.htm

    """

    linecolors = {'firebrick1': np.array([255, 48, 48])/255.,
                  'firebrick': np.array([178, 34, 34])/255.,
                  'orangered': np.array([255, 69, 0])/255.,
                  'orangered4': np.array([139, 37, 0])/255.,
                  'darkgoldenrod1': np.array([255, 185, 15])/255.,
                  'darkgoldenrod': np.array([184, 134, 11])/255.,
                  'darkolivegreen1': np.array([202, 255, 112])/255.,
                  'darkolivegreen3': np.array([162, 205, 90])/255.,
                  'darkseagreen': np.array([143, 188, 143])/255.,
                  'darkseagreen4': np.array([105, 139, 105])/255.,
                  'mediumpurple1': np.array([171, 130, 255])/255.,
                  'mediumpurple4': np.array([93, 71, 139])/255.,
                  'dodgerblue': np.array([30, 144, 255])/255.,
                  'mediumblue': np.array([0, 0, 205])/255.,
                  'limegreen': np.array([50, 205, 50])/255., 
                  'lightsteelblue3': np.array([162, 181, 205])/255., # more grey looking
                  'lightsteelblue4': np.array([110, 123, 139])/255., # more grey looking
                  'steelblue3': np.array([79, 148, 205])/255.,  # more blue looking
                  'steelblue4': np.array([54, 100, 139])/255.,
                  'pink2': np.array([238, 169, 184])/255.,
                  'pink4': np.array([139, 99, 108])/255.,
                  'chocolate1': np.array([255, 127, 36])/255.,
                  'chocolate4': np.array([139, 69, 19])/255.,
                  'barnred': np.array([124,10,2])/255.,
                  'yelloworange': np.array([255, 155, 0])/255.,
                  'darkyellow': np.array([242, 228, 40])/255.,
                  'magenta': np.array([112, 0, 102])/255.,
                  'violet': np.array([85, 0, 140])/255.,
                  'deepskyblue': np.array([0,191,255])/255.,
                  'midnightblue': np.array([25,25,112])/255.,
                  'skyblue': np.array([135, 206, 250])/255.,
                  'red1': np.array([118, 1, 1])/255.,    # dark to light: these reds are kind of pinky
                  'red2': np.array([171, 31, 31])/255.,
                  'red3': np.array([205, 72, 72])/255.,
                  'red4': np.array([228, 121, 121])/255.,
                  'red5': np.array([252, 191, 191])/255.,
                  'orange1': np.array([145, 73, 0])/255.,  # dark to light: these oranges are browny
                  'orange2': np.array([185, 93, 1])/255.,
                  'orange3': np.array([219, 137, 35])/255.,
                  'orange4': np.array([246, 144, 41])/255.,
                  'orange5': np.array([250, 193, 135])/255.,
                  'warm1': np.array([150, 1, 1])/255.,
                  'warm2': np.array([215, 50, 1])/255., 
                  'warm3': np.array([240, 119, 1])/255.,
                  'warm4': np.array([255, 163, 1])/255.,
                  'warm5': np.array([255, 200, 1])/255.,
                  'turq': np.array([0, 153, 153])/255.,
                  'niceblue': np.array([0,102,204])/255.,
                  'niceblue2': np.array([0,102,175])/255.,
                  'paperblue': np.array([0,56,226])/255.}

    return linecolors

def get_colordict(project=None):
    """ get_colordict(project=None):
                            returns a dictionary of line colors
                            by simulation key for the prescribed sea-ice
                            AGCM simulations.

                            project: not implemented yet, but there as a
                                     placeholder for when I add new project
                                     simulations
                            
                            Current keys:
                            

                            """

    #pairs = con.get_simpairsdict()
    ##     pairs.keys()
    ## Out[33]: 
    ## ['NSIDC',
    ##  'R4',
    ##  'CANnothk',
    ##  'R1',
    ##  'R2',
    ##  'R5',
    ##  'RCPa',
    ##  'R3',
    ##  'HAD',
    ##  'ENS',
    ##  'CAN',
    ##  'CANnosst',
    ##  'R4ct']

    
    colordict = {'R1': get_linecolor('warm1'), #'firebrick'),
             'R4': get_linecolor('warm2'), #'firebrick1'),
             'R3': get_linecolor('warm3'), #'yelloworange'),#'chocolate1'),
             'R5': get_linecolor('warm4'), #'darkyellow') #'skyblue'), #yelloworange'),
             'R2': get_linecolor('warm5'), #'steelblue3'), #'darkgoldenrod1'),
             'ENS': get_linecolor('magenta'),
             'CAN': get_linecolor('mediumblue'), #'mediumpurple1'), #darkyellow'),                 
             'HAD': get_linecolor('deepskyblue'),
             'NSIDC': get_linecolor('turq'),
             'R4ct': get_linecolor('darkseagreen4'),
             'CANnosst': get_linecolor('darkolivegreen3'),
             'CANnothk': get_linecolor('darkolivegreen1'),
             'RCPa': get_linecolor('orange1'),
             'E1': get_linecolor('mediumblue'),# same as CAN
             'E2': get_linecolor('mediumpurple1'),
             'E3': get_linecolor('mediumpurple4'),
             'E4': get_linecolor('violet'),
             'E5':  get_linecolor('lightsteelblue4'),
             'ENSE': get_linecolor('midnightblue'),
             'ESPR': get_linecolor('limegreen')}

    return colordict

def show_linecolors():

    linecolors = get_linecolorwheel()
    xx = np.array([0, 50])
    
    plt.figure(figsize=(5,10))
    lckeys=linecolors.keys()
    lckeys.sort()
    for cii,ckey in enumerate(lckeys):
        #print cii
        #print ckey
        plt.plot(xx,np.array([cii, cii]),color=linecolors[ckey],linewidth=4)
        plt.text(1,cii-.5,ckey,fontsize=9)

    plt.ylim((-1,len(linecolors)))

    ## legend(('firebrick1','firebrick','orangered','orangered4',
    ##        'darkgoldenrod1','darkgoldenrod','lightsteelblue3','lightsteelblue4', ...
    ##     'darkolivegreen1','darkolivegreen3','darkseagreen','darkseagreen4', ...
    ##     'mediumpurple1','mediumpurple4','steelblue3','steelblue4', ...
    ##     'dodgerblue','mediumblue','limegreen',...
    ##     'location','eastoutside')



        
# this means, if we are running this module as main script
#    (not importing from another), register and show
if __name__ == "__main__":
    register_cccmacms()
    show_cccmacms()
    show_linecolors()
    print 'http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps'
    print 'http://matplotlib.org/1.2.1/examples/pylab_examples/show_colormaps.html'
else:
    register_cccmacms()
