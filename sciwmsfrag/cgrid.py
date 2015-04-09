'''
Created on Apr 9, 2015

@author: ayan
'''
import numpy as np
from matplotlib.cm import get_cmap


def subset(latmin, lonmin, latmax, lonmax, lat, lon):
    #t1 = timeobj.time()
    latbool = (lat <= latmax+.18) & (lat >= latmin-.18)
    lonbool = (lon <= lonmax+.18) & (lon >= lonmin-.18)
    index = np.asarray(np.where(latbool & lonbool)).squeeze()
        #((lat <= latmax) == (lat >= latmin)) ==
        #((lon <= lonmax) == (lon >= lonmin),))).squeeze()
    #if (lonmax > 0) & (lonmin < 0):
    #    lon[lon > lonmax+30] = np.nan # would prefer to be subsetting the smallest area possible isntead of just hacking the rendering...
    #    lon[lon < lonmin-30] = np.nan
    if index.shape[1] > 0:
        ind = np.asarray(range(np.min(np.min(index[0])),np.max(np.max(index[0]))+1))
        jnd = np.asarray(range(np.min(np.min(index[1])),np.max(np.max(index[1]))+1))
        lat = lat[ind[0]:ind[-1], jnd[0]:jnd[-1]]
        lon = lon[ind[0]:ind[-1], jnd[0]:jnd[-1]]
    else:
        index = None
        lat = np.asarray([[],[]])
        lon = np.asarray([[],[]])
    #print str(timeobj.time()-t1) + " subset coords"
    return index, lat, lon


def plot(lon, lat, var1, var2, actions, ax, fig, **kwargs):
    #t1 = timeobj.time()
    aspect = kwargs.get('aspect', None)
    height = kwargs.get('height')
    width = kwargs.get('width')
    norm = kwargs.get('norm')
    blah = kwargs.get('cmap', 'jet')
    print('blah: {0}'.format(blah))
    cmap = get_cmap('jet')
    if "vectors" in actions:
        fig.set_figheight(height/80.0/aspect)
        fig.set_figwidth(width/80.0)
        vectors(lon, lat, var1, var2, mag, ax, norm, cmap, magnitude)
            
            
def vectors(lon, lat, var1, var2, mag, ax, norm, cmap, magnitude):
    if magnitude == "True":
        arrowsize = None
    elif magnitude == "False":
        arrowsize = 2.
    elif magnitude == "None":
        arrowsize = None
    else:
        arrowsize = float(magnitude)
    stride = 1
    ax.quiver(lon[::stride,::stride], lat[::stride,::stride], var1.squeeze()[::stride,::stride], var2.squeeze()[::stride,::stride], mag.squeeze()[::stride,::stride],
                pivot='mid',
                #units='uv', #xy
                cmap=cmap,
                norm=norm,
                minlength=.5,
                scale=arrowsize,
                scale_units='inches',
                angles='uv',
                )