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
VERTICAL_SLICE = -1  # get the last vertical slice

os.environ['TCL_LIBRARY'] = 'C:/Python279/tcl/tcl8.5'
os.environ['TK_LIBRARY'] = 'C:/Python279/tcl/tk8.5'


def display_shape(var_array, display_text):
    array_shape = var_array.shape
    message = '{display_text} shape: {shape}'.format(display_text=display_text,
                                                     shape=array_shape
                                                     )
    print(message)


if __name__ == '__main__':

    coawst = nc4.Dataset(SGRID_URL)
    sgc = SGrid().from_nc_dataset(coawst)
    grid_center_lon = sgc.grid_cell_center_lon
    grid_center_lat = sgc.grid_cell_center_lat
    u_raw = coawst.variables['u'][TIME_SLICE, VERTICAL_SLICE, :, :]
    display_shape(u_raw, 'u raw')
    v_raw = coawst.variables['v'][TIME_SLICE, VERTICAL_SLICE, :, :]
    display_shape(v_raw, 'v raw')
    angle = coawst.variables['angle'][:]
    display_shape(angle, 'angle')
    rho_mask = coawst.variables['mask_rho'][:]
    display_shape(rho_mask, 'rho mask')
    print(rho_mask)
    face_padding = sgc.face_padding
    print(face_padding)
    u_trim = u_raw[1:-1, :]
    v_trim = v_raw[:, 1:-1]
    display_shape(u_trim, 'u trim')
    display_shape(v_trim, 'v trim')
    angle_trim = angle[1:-1, 1:-1]
    rho_mask_trim = rho_mask[1:-1, 1:-1]
    display_shape(angle_trim, 'angle trim')
    display_shape(rho_mask_trim, 'rho mask trim')