# Import necessary libraries
import numpy as np
from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from draw_turbine import draw_turbine
yd = {'rotation': 0, 'ha': 'right', 'va': 'center'}

'''
    Functions to plot post-processed data from farm2d.
'''

# Contour plot (for unstructured and structured data)
def plot_xyslice(X,Y,VAR,varname='$U$ [m/s]',varlim='None',lvls=50,filename='None'):
    if varlim == 'None':
        varmin = np.min(VAR); varmax = np.max(VAR)
    else:
        varmin = varlim[0]; varmax = varlim[1];
    lvls = np.linspace(varmin, varmax, lvls)
    fig = plt.figure(figsize=(8,4))
    if len(np.shape(X)) == 1:
        # Unstructured data in vector form
        p = plt.tricontourf(X, Y, VAR, lvls, cmap=cm.turbo)
    elif len(np.shape(X)) == 2:
        # Structured data in matrix form
        p = plt.contourf(X, Y, VAR, lvls, cmap=cm.turbo)
    plt.xlabel('$x$ [m]')
    plt.ylabel('$y$ [m]',**yd)
    plt.axis('scaled')
    divider = make_axes_locatable(fig.gca())
    cax = divider.append_axes("right", size="4%", pad=0.1)
    comcolRANS = plt.colorbar(p, cax=cax)
    cax.set_ylabel(varname,rotation=0,va='center',ha='left')
    comcolRANS.set_ticks(np.linspace(varmin,varmax,5))
    if filename is not 'None':
        fig.savefig(filename,bbox_inches='tight')
    return fig


def get_xline(sdata, var, x):
    """
        For structured data, get line with constant x.
        Finds the closest x-line. One can interpret this as an "step"-interpolation
        , i.e. where the variable is constant within each cell.
    """
    X = sdata['X']
    Y = sdata['Y']
    U = sdata[var]
    # Find the column index in X
    x_idx = np.argmin(np.abs(X[0, :] - x))
    # Extract the y-coordinates and U values at this column index
    y_line = Y[:, x_idx]
    u_line = U[:, x_idx]
    return y_line, u_line

def get_yline(sdata, var, y):
    X = sdata['X']
    Y = sdata['Y']
    U = sdata[var]
    y_idx = np.argmin(np.abs(Y[:, 0] - y))
    x_line = X[y_idx, :]
    u_line = U[y_idx, :]
    return u_line

def get_ADflow(sdata,var,yADpos,D):
    '''
        Calculate disk-averaged flow.
        This simple calculation is only valid for uniform cells.
    '''
    # Assuming X, Y, U are your input 2D arrays, and ymin, ymax are your scalar values
    ymin = yADpos*D - D/2  # Set your desired ymin value
    ymax = yADpos*D + D/2  # Set your desired ymax value

    Y = sdata['Y']
    U = sdata[var]
    yminidx = np.argmin(np.abs(Y[:, 0] - ymin))
    ymaxidx = np.argmin(np.abs(Y[:, 0] - ymax))
    Ydisk = Y[yminidx:ymaxidx+1,:]
    Udisk = U[yminidx:ymaxidx+1,:]
    Udisk = np.sum(Udisk,axis=0)/np.shape(Udisk)[0]
    return Udisk

def plot_longwake(Xn,Yn,Un,xval,xmin,xmax):
    '''
        Coordinates should be normalized by D. Also, turbine should be at (0,0).
        Un = 1 - U/Uinf.
        xval is the downstream positions (in D), where lines are extracted.
    '''
    fig = plt.figure(figsize=(12,3))
    ax = plt.gca()
    data = {'X': Xn, 'Y': Yn, 'Un': Un}
    for x_val in xval:
        y_line, u_line = get_xline(data, 'Un', x_val)
        plt.plot(u_line + x_val, y_line, 'b-')

    secax = ax.secondary_xaxis('top')
    secax.set_xlabel("$x/D$")
    secax.set_xticks(np.append(np.insert(xval, 0, xmin),xmax))

    def half_array(x):
        new_array = []
        one_array = []
        for i in range(len(x)):
            new_array.append(x[i])
            new_array.append(x[i]+0.5)
            one_array.append(1.0)
            one_array.append(0.5)
        return new_array, one_array

    new_array, one_array = half_array(xval)
    ax.set_xticks(new_array)
    ax.set_xticklabels(one_array)

    ax.set_xlabel(r'$U/U_\infty$')
    ax.set_ylabel('$y/D$',**yd)

    draw_turbine(D=1)

    plt.ylim([-1.5,1.5])
    plt.xlim([xmin-0.5,xmax+0.5])