# Import necessary libraries
import numpy as np
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


# Load whole iteration history OpenFOAM AD output file.
def load_ad_history(filename):
    # Get headers
    with open(filename, 'r') as file:
        #file.readline()  # Skip the first line
        #headers = file.readline().strip().split()  # Read the second line for headers
        file.readline()  # Skip the first line
        headers = file.readline().strip().split()
        headers = [h for h in headers if h != '#']  # Remove all '#' symbols

    # Read actually data
    data = np.genfromtxt(filename,skip_header=2)

    # Create a dictionary to store all iterations of each column
    disk = {header: data[:, i] for i, header in enumerate(headers)}

    return disk


# Load all ADs.
def load_all_ad(turbine_file='turbines.dat',post_folder='postProcessing',twod=True,rho=1.225):
    turbines = np.genfromtxt(turbine_file,skip_header=1)
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
        diski.update(load_ad('%s/disk%d/0/AD_calaf.dat'%(post_folder,wt_nr[i])))

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

def extract_number_of_cores(filename):
    nprocs = None
    with open(filename, 'r') as file:
        for line in file:
            if "nProcs :" in line:
                parts = line.split()
                try:
                    nprocs = float(parts[2])
                except (IndexError, ValueError):
                    pass
    return int(nprocs)

def extract_number_of_iterations(filename):
    iterations = None
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('Time = '):
                parts = line.split()
                try:
                    iterations = float(parts[2])
                except (IndexError, ValueError):
                    pass
    return int(iterations)

def extract_info(filename):
    walltime = extract_last_execution_time(filename)
    cores = extract_number_of_cores(filename)
    iterations = extract_number_of_iterations(filename)
    cputime = walltime*cores
    res = {'walltime': walltime,
           'cores': cores,
           'iterations': iterations,
           'cputime': cputime,
           }
    return res


def extract_number_of_cells(filename):
    """
        Load number of cells from log.blockMesh file.
    """
    ncells = None
    with open(filename, 'r') as file:
        for line in file:
            if "nCells" in line:
                parts = line.split()
                try:
                    ncells = float(parts[1])
                except (IndexError, ValueError):
                    pass
    return int(ncells)


