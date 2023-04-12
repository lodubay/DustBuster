"""
Create an interactive 3D plot of the Lallement et al. (2022) differential
extinction map.
"""

from astropy.io import fits
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from coordinate_transform_functions import radec_to_data

# plot range in pc
# [(xmin, xmax), (ymin, ymax), (zmin, zmax)]
plot_range = [(-3005, 2995), (-3005, 2995), (-405, 395)]
# location of Sun (origin) in data grid units
sun_pos = (300.5, 300.5, 40.5)
# parsecs per grid step
gridstep_pc = 10
# number of indices to skip between samples in (x, y, z) axes
sampling = 4
# coordinates of Cygnus OB2
CygOB2_coords = (np.array(radec_to_data(308.3, 41.32, 1570)) - np.array(sun_pos)) * gridstep_pc
print(CygOB2_coords)

# import data
hdul = fits.open('data/cube_ext.fits')
data = hdul[0].data
# invert data axes
data = np.swapaxes(data, 0, 2)
# downsample
data = data[::sampling, ::sampling, ::sampling]
gridstep_pc *= sampling
sun_pos = (sun_pos[0]/sampling, sun_pos[1]/sampling, sun_pos[2]/sampling)
print(data.shape)

# fix issues with user-defined plot range
for i in range(3):
    # adjust to proper spacing
    # plot_range[i] =  (round(plot_range[i][0] - 5, -1) + 5,
    #                   round(plot_range[i][1] - 5, -1) + 5)
    # maximum plot range in given axis
    plot_range[i] = (max(plot_range[i][0], 
                         -sun_pos[i] * gridstep_pc),
                     min(plot_range[i][1], 
                         (data.shape[i] - 1 - sun_pos[i]) * gridstep_pc))
print(plot_range)

# convert plot range to data cube slice
bounds = []
for i in range(3):
    # ensure plot range does not exceed data limits
    xmin = int(plot_range[i][0] / gridstep_pc + sun_pos[i])
    xmax = int(plot_range[i][1] / gridstep_pc + sun_pos[i]) + 1
    bounds.append(slice(max(xmin, 0), 
                        min(xmax, data.shape[i])))
print(bounds)
# subset of data cube
subset = data[tuple(bounds)]

# make meshgrid
xvals = np.arange(plot_range[0][0], plot_range[0][1]+gridstep_pc, gridstep_pc)
yvals = np.arange(plot_range[1][0], plot_range[1][1]+gridstep_pc, gridstep_pc)
zvals = np.arange(plot_range[2][0], plot_range[2][1]+gridstep_pc, gridstep_pc)
X, Y, Z = np.meshgrid(xvals, yvals, zvals, indexing='ij')
print(X.shape)

# volume plot of extinction
fig = go.Figure(
    data=[
        go.Volume(
            x=X.flatten(),
            y=Y.flatten(),
            z=Z.flatten(),
            value=subset.flatten(),
            isomin=1e-3,
        #     isomax=0.8,
            opacity=0.1, # needs to be small to see through all surfaces
            surface_count=20, # needs to be a large number for good volume rendering
            colorbar=dict(
                title='Extinction density [nanomag/pc]'
            )
        ),
        # plot location of Cygnus OB2
        go.Scatter3d(
            x=CygOB2_coords[0:1],
            y=CygOB2_coords[1:2],
            z=CygOB2_coords[2:],
            mode='markers'
        )
    ]
)

fig.update_layout(
    scene=dict(
        xaxis_title='x [pc]',
        yaxis_title='y [pc]',
        zaxis_title='z [pc]'
    ),
    scene_aspectmode='data',
    title=dict(
        text='Volumetric Extinction Density Map (Lallement et al. 2022)'
    )
)

fig.write_html('3d_dust_map.html')
