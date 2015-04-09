'''
Created on Apr 9, 2015

@author: ayan
'''
import os
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize
from matplotlib.backends.backend_agg import FigureCanvasAgg
from mpl_toolkits.basemap import Basemap
from pysgrid import SGrid
from pysgrid.processing_2d import vector_sum, rotate_vectors, avg_to_cell_center
from render_map import SGRID_URL
from sciwmsfrag.cgrid import subset


CURRENT_DIR = os.path.dirname(__name__)
CACHE_FILE = os.path.join(CURRENT_DIR, 'coawst_cache.nc')


if __name__ == '__main__':
    
    lonmin = -101.74989
    lonmax = -53.253086
    latmin = 11.888468
    latmax = 48.463924
    variables = ['u', 'v']
    ds = nc4.Dataset(SGRID_URL)
    sg = SGrid()
    ds_sgrid = sg.from_nc_dataset(ds)
    # ds_sgrid.save_as_netcdf(CACHE_FILE)
    cache_ds = nc4.Dataset(CACHE_FILE)
    lon_centers = cache_ds.variables['grid_center_lon'][:]
    lat_centers = cache_ds.variables['grid_center_lat'][:]
    subsetting = subset(latmin, lonmin, latmax, lonmax, lat_centers, lon_centers)
    index, lat, lon = subsetting
    times = cache_ds.variables['time'][:]
    time_idx = -1
    layer_idx = 0
    fig = Figure(dpi=80, facecolor='blue', edgecolor='green')
    fig.set_alpha(0)
    m = Basemap(llcrnrlon=lonmin, 
                llcrnrlat=latmin,
                urcrnrlon=lonmax, 
                urcrnrlat=latmax,
                resolution=None,
                lat_ts = 0.0,
                suppress_ticks=True,
                projection='merc'
                )
    m.ax = fig.add_axes([0, 0, 1, 1], xticks=[], yticks=[])
    tl = np.s_[time_idx, layer_idx]
    u = ds.variables['u'][tl + ds_sgrid.u_slice[2:]]
    print(u.shape)
    v = ds.variables['v'][tl + ds_sgrid.v_slice[2:]]
    print(v.shape)
    lonc = lon_centers[ds_sgrid.lon_rho_slice]
    latc = lon_centers[ds_sgrid.lat_rho_slice]
    u_avg = avg_to_cell_center(u, 1)  # x-directed velocity
    v_avg = avg_to_cell_center(v, 0)  # y-directed velocity
    angles = ds.variables['angle'][ds_sgrid.angle_slice]
    u_rot, v_rot = rotate_vectors(u_avg, v_avg, angles)
    uv_mag = vector_sum(u_rot, v_rot)
    stride = 1
    cmap = get_cmap('jet')
    norm = Normalize()
    height = 800
    width = 1200
    fig.set_figheight(height/80.0/m.aspect)
    fig.set_figwidth(width/80.0)
    m.ax.quiver(lonc[::stride, ::stride],
                latc[::stride, ::stride],
                u_rot.squeeze()[::stride, ::stride],
                v_rot.squeeze()[::stride, ::stride],
                uv_mag.squeeze()[::stride, ::stride],
                pivot='mid',
                cmap=cmap,
                norm=norm,
                minlength=0.5,
                scale=None,
                scale_units='inches',
                angles='uv'
                )
    lonmax, latmax = m(lonmax, latmax)
    lonmin, latmin = m(lonmin, latmin)
    m.ax.set_xlim(lonmin, lonmax)
    m.ax.set_ylim(latmin, latmax)
    m.ax.set_frame_on(False)
    m.ax.set_clip_on(False)
    m.ax.set_position([0, 0, 1, 1])
    canvas = FigureCanvasAgg(fig)
    canvas.print_png(r'C:\Users\ayan\Desktop\tmp\test.png')
    
    