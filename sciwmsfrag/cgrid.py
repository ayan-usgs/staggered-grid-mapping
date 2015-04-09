'''
Created on Apr 9, 2015

@author: ayan
'''
import numpy as np


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