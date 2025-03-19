from farm2d import farm2d
from sys_commands import *
from farm2d_post import *
from farm2d_plot import *
from farm2d_pyvista import *

'''
    Horns Rev 1 simulation.
    farm2d: https://github.com/mchba/farm2d
'''

clean = 1
files = 1
run   = 1
post  = 1


####################### Initialize ###############################
# From https://gitlab.windenergy.dtu.dk/TOPFARM/PyWake/-/blob/master/py_wake/examples/data/hornsrev1.py?ref_type=heads
wt_x = [423974, 424042, 424111, 424179, 424247, 424315, 424384, 424452, 424534,
        424602, 424671, 424739, 424807, 424875, 424944, 425012, 425094, 425162,
        425231, 425299, 425367, 425435, 425504, 425572, 425654, 425722, 425791,
        425859, 425927, 425995, 426064, 426132, 426214, 426282, 426351, 426419,
        426487, 426555, 426624, 426692, 426774, 426842, 426911, 426979, 427047,
        427115, 427184, 427252, 427334, 427402, 427471, 427539, 427607, 427675,
        427744, 427812, 427894, 427962, 428031, 428099, 428167, 428235, 428304,
        428372, 428454, 428522, 428591, 428659, 428727, 428795, 428864, 428932,
        429014, 429082, 429151, 429219, 429287, 429355, 429424, 429492]
wt_y = [6151447, 6150891, 6150335, 6149779, 6149224, 6148668, 6148112, 6147556,
        6151447, 6150891, 6150335, 6149779, 6149224, 6148668, 6148112, 6147556,
        6151447, 6150891, 6150335, 6149779, 6149224, 6148668, 6148112, 6147556,
        6151447, 6150891, 6150335, 6149779, 6149224, 6148668, 6148112, 6147556,
        6151447, 6150891, 6150335, 6149779, 6149224, 6148668, 6148112, 6147556,
        6151447, 6150891, 6150335, 6149779, 6149224, 6148668, 6148112, 6147556,
        6151447, 6150891, 6150335, 6149779, 6149224, 6148668, 6148112, 6147556,
        6151447, 6150891, 6150335, 6149779, 6149224, 6148668, 6148112, 6147556,
        6151447, 6150891, 6150335, 6149779, 6149224, 6148668, 6148112, 6147556,
        6151447, 6150891, 6150335, 6149779, 6149224, 6148668, 6148112, 6147556]
D = 80

# Translate layout
wt_x = wt_x - np.min(wt_x) + 2*D
wt_y = wt_y - np.min(wt_y) + 2*D

# Quick plot of the layout
fig = plt.figure()
plt.plot(wt_x/D,wt_y/D,'bo')
plt.axis('scaled')
plt.xlabel('$x/D$')
plt.ylabel('$y/D$')
plt.xlim(left=0)
plt.ylim(bottom=0)

# Create turbines
def V80(x,y):
    return {'x': x,
            'y': y,
            'D': D,
            'zh': 70,
            'CT': 0.75
            # No CP needed, because 1D mom controller is used.
            }

wts = []
for i in range(len(wt_x)):
    wts.append(V80(wt_x[i],wt_y[i]))

par = {'D': D,
  'lx': 75*D,
  'ly': 55*D,
  'dx': D/8,
  'Lxw': 50*D,
  'Lxe': 50*D,
  'Lyn': 50*D,
  'Lys': 50*D,
  'dxo': 4*D,
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
    plot_xyslice(sdata['X'], sdata['Y'], sdata['U'], varname='$U$ [m/s]',lvls=40,filename='hornsrev1_u_contour.png')

