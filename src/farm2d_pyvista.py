# Import necessary libraries
import numpy as np
import xarray as xr
import pyvista as pv

'''
    Functions to load flow field using PyVista.
'''

# Read OpenFOAM 2D xy-flowfield with pyvista.
def load_flowfield(foamfile):
    reader = pv.get_reader(foamfile)
    reader.set_active_time_value(reader.time_values[-1])  # use last time point
    mesh = reader.read()
    internalMesh = mesh['internalMesh']
    boundaries = mesh['boundary']

    # Grid
    grid = internalMesh.points
    x = grid[:,0]
    y = grid[:,1]
    outdata = {'x': x, 'y': y}

    # Data
    data = internalMesh.point_data
    fields = list(data.keys())
    for k in fields:
        if k == 'U':
            outdata['U'] = data['U'][:,0]
            outdata['V'] = data['U'][:,1]
        if k == 'p':
            outdata['P'] = data['p']
        if k == 'k':
            outdata['K'] = data['k']
        if k == 'epsilon':
            outdata['eps'] = data['epsilon']
        if k == 'turbulenceProperties:R':
            outdata['uu'] = data['turbulenceProperties:R'][:,0] # probably correct
            outdata['vv'] = data['turbulenceProperties:R'][:,1] # should probably be vv
            outdata['ww'] = data['turbulenceProperties:R'][:,2] # should probably be ww
            outdata['uv'] = data['turbulenceProperties:R'][:,3] # should probably be uv
            outdata['uw'] = data['turbulenceProperties:R'][:,4] # should proably be uw
            outdata['vw'] = data['turbulenceProperties:R'][:,5] # should probably be vw

    return outdata


# Read a structured rectilinear OpenFOAM flowfield
def load_sflowfield(foamfile):
    # Load unstructured data
    udata = load_flowfield(foamfile)
    x = udata['x']
    y = udata['y']

    # Create 2D arrays (can do this if the grid is structured and rectilinear)
    x_unique = np.unique(x)
    y_unique = np.unique(y)

    # Create meshgrid from unique x and y
    X, Y = np.meshgrid(x_unique, y_unique)

    # Initialize an empty 2D array for U
    sdata = {'x': x_unique, 'y': y_unique, # 1D arrays
            'X': X, 'Y': Y   # 1D arrays
            }

    # Slow version
    # for k in udata.keys():
    #     if k == 'x' or k == 'y':
    #         continue
    #     var = udata[k]
    #     var_2D = np.zeros((len(y_unique), len(x_unique)))
    #     # Fill the var_2D array by mapping var's 1D values to the correct (x, y) positions
    #     for i in range(len(x)):
    #         x_idx = np.where(x_unique == x[i])[0][0]  # Find the x index in the unique values
    #         y_idx = np.where(y_unique == y[i])[0][0]  # Find the y index in the unique values
    #         var_2D[y_idx, x_idx] = var[i]  # Place the var value at the correct position
    #         sdata[k] = var_2D

    # Faster version
    for k in udata.keys():
        if k == 'x' or k == 'y':
            continue
        sdata[k] = np.zeros((len(y_unique), len(x_unique)))
    # Fill the var_2D array by mapping var's 1D values to the correct (x, y) positions
    for i in range(len(x)):
        x_idx = np.where(x_unique == x[i])[0][0]  # Find the x index in the unique values
        y_idx = np.where(y_unique == y[i])[0][0]  # Find the y index in the unique values
        for k in udata.keys():
            if k == 'x' or k == 'y':
                continue
            sdata[k][y_idx, x_idx] = udata[k][i]  # Place the var value at the correct position

    return sdata

# Read flowfield into a netCDF file
def load_nflowfield(foamfile):
    # Load the structured data
    sdata = load_sflowfield(foamfile)

    dv = {}
    for k in sdata.keys():
        if k=='x' or k=='y' or k=='X' or k=='Y':
            continue
        dv[k] = (["x", "y"], sdata[k].T)

    # Convert to netCDF format
    data = xr.Dataset(
    data_vars=dv,
    coords=dict(
        x=("x", sdata['x']),
        y=("y", sdata['y'])
    ),
    attrs=dict(case='V80')
    )
    return data