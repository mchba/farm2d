from farm2d import farm2d
from sys_commands import *
from farm2d_post import *
from farm2d_plot import *
from farm2d_pyvista import *

'''
    Single V80 turbine simulation.
    farm2d: https://github.com/mchba/farm2d
'''

clean = 1
files = 1
run   = 1
post  = 1


####################### Initialize ###############################
D = 80
def V80(x,y):
    return {'x': x,
            'y': y,
            'D': D,
            'zh': 70,
            'CT': 0.75
            # No CP needed, because 1D mom controller is used.
            }

wts = [V80(3*D,2*D),
#       V80(10*D,2*D)
   ]

par = {'D': D,
  'lx': 16*D,
  'ly': 4*D,
  'dx': D/8,
  'Lxw': 40*D,
  'Lxe': 40*D,
  'Lyn': 25*D,
  'Lys': 25*D,
  'dxo': 2*D,
  'Uinf': 8.0,
  'Iinf': 0.054,
  'zref': 70,
  'wts': wts
  }

F2D = farm2d()
F2D.setup(par)


####################### Clean folder from old run-files ################
if clean == 1:
    system_com('./Allclean')


####################### Setup OpenFOAM files ################
if files == 1:
    F2D.make_grid('system/blockMeshDict')
    F2D.plot_grid(save_fig=True)
    F2D.make_inflow('0/ASLparameters')
    F2D.make_topo('system/topoSetDict')
    F2D.make_fvoptions('constant/fvOptions')
    F2D.make_turbulence('constant/turbulenceProperties')


############# Run simulation #########################
if run == 1:
    system_com('./Allrun')


############# Post-process #########################
if post == 1:
    foamfile = 'test.foam'
    sdata = load_sflowfield(foamfile)
    plot_xyslice(sdata['X'], sdata['Y'], sdata['U'], varname='$U$ [m/s]',lvls=20,filename='farm2d_u_contour.png')

