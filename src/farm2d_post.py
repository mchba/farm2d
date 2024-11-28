# Import necessary libraries
import numpy as np
import pyvista as pv
import xarray as xr

'''
    Functions to load and post-process data from farm2d.
'''

# Load data from OpenFOAM AD output file.
def load_ad(filename):
    # Get headers
    with open(filename, 'r') as file:
        #file.readline()  # Skip the first line
        #headers = file.readline().strip().split()  # Read the second line for headers
        file.readline()  # Skip the first line
        headers = file.readline().strip().split()
        headers = [h for h in headers if h != '#']  # Remove all '#' symbols
        
    # Read actually data
    data = np.genfromtxt(filename,skip_header=2)
    
    # Create a dictionary to store final value of each column
    disk = {header: data[-1, i] for i, header in enumerate(headers)}
    
    return disk


# Load all ADs.
def load_all_ad(twod=True,rho=1.225):
    turbines = np.genfromtxt('turbines.dat',skip_header=1)
    if len(np.shape(turbines)) == 1:
        # Special case: If there is only one turbine, convert 1D array to 2D array.
        turbines = np.reshape(turbines, (-1,np.shape(turbines)[0]))
    wt_nr = turbines[:,0]
    wt_x = turbines[:,1]
    wt_y = turbines[:,2]
    wt_zh = turbines[:,3]
    wt_D = turbines[:,4]
    
    # List of AD dictionaries
    turbines_out = []
    for i in range(len(wt_x)):
        # Basic setup
        diski = {'disknr': int(wt_nr[i]),
                'x':  wt_x[i],
                'y':  wt_y[i],
                'zh': wt_zh[i],
                'D':  wt_D[i],
                 }
        # Combine with OpenFOAM output
        diski.update(load_ad('postProcessing/disk%d/0/AD_calaf.dat'%(wt_nr[i])))
        
        # Correct T = 0.5*rho*D*U**2*CT to T = 0.5*rho*A*U**2*CT
        if twod == True:
            diski['T'] = diski['T']*np.pi*diski['D']/4
            diski['P'] = diski['P']*np.pi*diski['D']/4
        
        # simpleFOAM assumes rho=1, but can multiply with rho in post-process to get real power.
        diski['T'] = diski['T']*rho
        diski['P'] = diski['P']*rho
        
        # Save to array
        turbines_out.append(diski)
    return turbines_out
    
    
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
            outdata['uu'] = data['turbulenceProperties:R'][:,0] 
            outdata['uv'] = data['turbulenceProperties:R'][:,1] 
            outdata['uw'] = data['turbulenceProperties:R'][:,2] 
            outdata['vv'] = data['turbulenceProperties:R'][:,3] 
            outdata['vw'] = data['turbulenceProperties:R'][:,4] 
            outdata['ww'] = data['turbulenceProperties:R'][:,5] 
            
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
    

# Load line data output by OpenFOAM
def load_line(filename):
    '''
        The below code is only valid for lines with (U p k nut) data.
        Should be generalized.
        However, it is easier to just load the whole structured flow field
        instead of using the OpenFOAM linesample approach.
    '''
    line = np.genfromtxt(filename)
    line = {'y': line[:,0],
            'K': line[:,1],
            'nut': line[:,2],
            'P': line[:,3],
            'U': line[:,4],
            'V': line[:,5],
            'W': line[:,6],
            }
    return line
    

def extract_last_execution_time(filename):
    """
        Load computational time from log.simpleFoam file.
    """
    execution_time = None
    with open(filename, 'r') as file:
        for line in file:
            if "ExecutionTime =" in line:
                parts = line.split()
                try:
                    execution_time = float(parts[2])
                except (IndexError, ValueError):
                    pass
    return execution_time
    
    
    
    