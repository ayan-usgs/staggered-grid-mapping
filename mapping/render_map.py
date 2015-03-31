'''
Created on Mar 25, 2015

@author: ayan
'''
import os
import netCDF4 as nc4
import numpy as np
from pysgrid.sgrid import SGrid


ADRIATIC_CGRID = 'http://geoport.whoi.edu/thredds/dodsC/examples/bora_feb.nc'
SGRID_URL = 'http://geoport.whoi.edu/thredds/dodsC/coawst_4/use/fmrc/coawst_4_use_best.ncd'
TIME_SLICE = -1  # get the last time slice

os.environ['TCL_LIBRARY'] = 'C:/Python279/tcl/tcl8.5'
os.environ['TK_LIBRARY'] = 'C:/Python279/tcl/tk8.5'



if __name__ == '__main__':

    coawst = nc4.Dataset(SGRID_URL)
    sgc = SGrid().from_nc_dataset(coawst)
    grid_center_lon = sgc.grid_cell_center_lon
    grid_center_lat = sgc.grid_cell_center_lat
    u_raw = coawst.variables['u'][:]
    v_raw = coawst.variables['v'][:]
    
