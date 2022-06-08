### setting for matplotlib
### Feb. 2022 / Masakiyo Kitazawa

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pylab import *
from matplotlib.colors import ListedColormap
import colorsys
from cycler import cycler

plt.rcParams.update({
    'figure.figsize' : [8,6] ,
    'font.size' : 20 ,
    'xtick.labelsize' : 18 ,
    'ytick.labelsize' : 18 ,
    'axes.labelsize' : 23 ,
    'axes.titlesize' : 23 ,

    'mathtext.fontset' : 'cm' ,

    'axes.grid' : True ,
    'grid.linestyle' : 'dotted' ,

    'axes.labelpad' : 10 ,
    'savefig.bbox' : 'tight' ,

    'lines.linewidth' : 3 ,
    'errorbar.capsize' : 6 ,

    #make linestyle 'dotted' denser
    'lines.dotted_pattern' : [1.0, 1.0]
})


#use darker color set
clist = []
for i in range(8):
    c = colorsys.rgb_to_hsv( *plt.get_cmap('Set1')(i)[:3] )
    c = colorsys.hsv_to_rgb( c[0] , 1.-(1.-c[1])*.2 , c[2]*.85 )
    # c = colorsys.rgb_to_hsv( *plt.get_cmap('tab10')(i)[:3] )
    # c = colorsys.hsv_to_rgb( c[0] , 1.-(1.-c[1])*.5 , c[2]*.95 )
    clist.append( c )
plt.rcParams['axes.prop_cycle'] = cycler( color=clist )

plt.gca().xaxis.get_major_formatter().set_useOffset(False)
plt.gca().yaxis.get_major_formatter().set_useOffset(False)

#for line styles
dashdot2 = ( 0 , (7,2,1,2,1,2) ) #dash-dot-dot

#same as xtics and ytics in gnuplot:
def xtics( tics ):
    plt.gca().xaxis.set_major_locator( ticker.MultipleLocator( tics ) )
    return
def ytics( tics ):
    plt.gca().yaxis.set_major_locator( ticker.MultipleLocator( tics ) )
    return

