# Import necessary libraries
import numpy as np
from sys_commands import system_com
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.backends.backend_pdf

def lineOut(id,s):
    line = ' '*id*4 + s + '\n'
    return line

class farm2d():
    '''
        Code to generate
          - blockMeshDict
          - topoSetDict
          - fvOptions
        for 2D OpenFOAM wind farm simulations.
    '''
    def __init__(self):
        '''
            Default parameters.
        '''

        # Turbine diameter used to normalize grid
        self.D = 80

        # Wake domain
        self.lx = 16*self.D
        self.ly = 4*self.D
        self.lz = 1
        self.dx = self.D/8

        # Outer domain
        self.dxo = 2*self.D  # Resolution at outer domain (used to stretch domain outwards)
        self.Lxw = 50*self.D    # West
        self.Lxe = 50*self.D    # East
        self.Lys = 50*self.D    # South
        self.Lyn = 50*self.D    # South

        # Turbines
        self.wts = []

        # Turbulence parameters
        self.Ce1  = 1.44
        self.Ce2  = 1.92
        self.Cmu  = 0.09
        self.sige = 1.3
        self.kap  = 0.41

        # Inflow (see Baungaard (2019, Eq.1.17-1.23) for summary of eqs)
        self.Uinf = 8.0
        self.Iinf = 0.054
        self.zref = 70
        self.Kinf = 1.5*(self.Iinf*self.Uinf)**2
        self.utau = self.Kinf**0.5*self.Cmu**0.25
        self.z0   = self.zref/(np.exp(self.kap*self.Uinf/self.utau) - 1)
        self.epsinf = self.utau**3/(self.kap*(self.zref + self.z0))

        # Source term
        self.kesource = True
        self.Sk = self.epsinf
        self.Seps = self.Ce2*self.epsinf**2/self.Kinf

    def setup(self,var):
        '''
            Method to load input settings.
        '''
        # Load variables from input dictionary.
        for k in var.keys():
            if hasattr(self, k):
                setattr(self, k, var[k])
            else:
                print(k + ' is not a valid input! Exiting...')

        # Update some variables
        self.Kinf = 1.5*(self.Iinf*self.Uinf)**2
        self.utau = self.Kinf**0.5*self.Cmu**0.25
        self.z0   = self.zref/(np.exp(self.kap*self.Uinf/self.utau) - 1)
        self.epsinf = self.utau**3/(self.kap*(self.zref + self.z0))
        self.Sk = self.epsinf
        self.Seps = self.Ce2*self.epsinf**2/self.Kinf


    def write_dict2file(self,s,filename,ob='unknown',ff=True):
        self.out = ''
        def l(id,s):
            self.out += lineOut(id,s)
        # Header
        l(0,"/*--------------------------------*- C++ -*----------------------------------*\\")
        l(0,"| =========                |                                                 |")
        l(0,"| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |")
        l(0,"|  \\    /   O peration     | Version:  v2206                                 |")
        l(0,"|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |")
        l(0,"|    \\/     M anipulation  |                                                 |")
        l(0,"\*---------------------------------------------------------------------------*/")
        if ff==True:
            l(0,"FoamFile")
            l(0,"{")
            l(1,"version         2.0;")
            l(1,"format          ascii;")
            l(1,"class           dictionary;")
            l(1,"object          %s;"%(ob))
            l(0,"}")
        l(0,"// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //")
        l(0,"")

        # Write to file
        with open(filename, "w") as file:
            file.write(self.out + s)


    def make_grid(self,filename):
        '''
            Create blockMeshDict
        '''

        # Block corners
        x1 = -self.Lxw
        x2 = 0.0
        x3 = self.lx
        x4 = self.lx+self.Lxe
        y1 = -self.Lys
        y2 = 0.0
        y3 = self.ly
        y4 = self.ly+self.Lyn
        z1 = 0.0
        z2 = self.lz

        # String to be made
        self.code_block = ''
        def l(id,s):
            self.code_block += lineOut(id,s)

        # Scale
        l(0,'scale 1;')
        l(0,'')

        # Vertices
        l(0,'vertices')
        l(0,'(')
        l(1,'(%.8f %.8f %.8f) // 0'%( x1,y1,z1))
        l(1,'(%.8f %.8f %.8f) // 1'%( x2,y1,z1))
        l(1,'(%.8f %.8f %.8f) // 2'%( x2,y2,z1))
        l(1,'(%.8f %.8f %.8f) // 3'%( x1,y2,z1))
        l(1,'(%.8f %.8f %.8f) // 4'%( x1,y1,z2))
        l(1,'(%.8f %.8f %.8f) // 5'%( x2,y1,z2))
        l(1,'(%.8f %.8f %.8f) // 6'%( x2,y2,z2))
        l(1,'(%.8f %.8f %.8f) // 7'%( x1,y2,z2))
        l(1,'(%.8f %.8f %.8f) // 8'%( x3,y1,z1))
        l(1,'(%.8f %.8f %.8f) // 9'%( x3,y2,z1))
        l(1,'(%.8f %.8f %.8f) // 10'%(x3,y1,z2))
        l(1,'(%.8f %.8f %.8f) // 11'%(x3,y2,z2))
        l(1,'(%.8f %.8f %.8f) // 12'%(x4,y1,z1))
        l(1,'(%.8f %.8f %.8f) // 13'%(x4,y2,z1))
        l(1,'(%.8f %.8f %.8f) // 14'%(x4,y1,z2))
        l(1,'(%.8f %.8f %.8f) // 15'%(x4,y2,z2))
        l(1,'(%.8f %.8f %.8f) // 16'%(x2,y3,z1))
        l(1,'(%.8f %.8f %.8f) // 17'%(x1,y3,z1))
        l(1,'(%.8f %.8f %.8f) // 18'%(x2,y3,z2))
        l(1,'(%.8f %.8f %.8f) // 19'%(x1,y3,z2))
        l(1,'(%.8f %.8f %.8f) // 20'%(x3,y3,z1))
        l(1,'(%.8f %.8f %.8f) // 21'%(x3,y3,z2))
        l(1,'(%.8f %.8f %.8f) // 22'%(x4,y3,z1))
        l(1,'(%.8f %.8f %.8f) // 23'%(x4,y3,z2))
        l(1,'(%.8f %.8f %.8f) // 24'%(x2,y4,z1))
        l(1,'(%.8f %.8f %.8f) // 25'%(x1,y4,z1))
        l(1,'(%.8f %.8f %.8f) // 26'%(x2,y4,z2))
        l(1,'(%.8f %.8f %.8f) // 27'%(x1,y4,z2))
        l(1,'(%.8f %.8f %.8f) // 28'%(x3,y4,z1))
        l(1,'(%.8f %.8f %.8f) // 29'%(x3,y4,z2))
        l(1,'(%.8f %.8f %.8f) // 30'%(x4,y4,z1))
        l(1,'(%.8f %.8f %.8f) // 31'%(x4,y4,z2))
        l(0,');')
        l(0,'')

        # Calculate stretching for all four directions
        def calc_stretch(dx1n,dx2n):
            # dx1n = dx1/L, dx2n = dx2/L
            r = (dx1n - 1)/(-1 + dx2n)
            N = int(np.ceil(np.log(dx2n/dx1n)/np.log(r) + 1))
            # Since we do a round-off with ceil (to get an integer N), actually the grid
            # is stretched slightly different than the stretching paramter r.
            rmod = np.exp(np.log(dx2n/dx1n)/(N-1))
            return N, rmod

        def stretch(x1,L,N,r):
            dx = np.zeros(N)
            dx[0] = 1  # We will scale it up later (we can not just use dx1, because N has been ceiled)
            for i in range(N-1):
                dx[i+1] = r*dx[i]

            x = np.zeros(N+1)
            x[0] = x1
            for i in range(N):
                x[i+1] = x[i]+dx[i]

            # Normalize x-array (make it go from 0 to 1)
            Lx = x[-1] - x[0]
            xnorm = (x-x[0])/Lx

            return xnorm*L+x[0]

        # Discretization (number of cells in the 9 blocks)
        nx1, rx1 = calc_stretch(self.dxo/self.Lxw, self.dx/self.Lxw)
        nx2 = self.lx/self.dx
        nx3, rx3 = calc_stretch(self.dx/self.Lxe, self.dxo/self.Lxe)
        ny1, ry1 = calc_stretch(self.dxo/self.Lys, self.dx/self.Lys)
        ny2 = self.ly/self.dx
        ny3, ry3 = calc_stretch(self.dx/self.Lyn, self.dxo/self.Lyn)
        nz1 = 1

        # Stretching (first to last cell in each block)
        Rw = self.dx/self.dxo
        Rn = self.dxo/self.dx
        Re = self.dxo/self.dx
        Rs = self.dx/self.dxo
        Rz = 1

        # Save rectilinear vertex coordinates (NOT cell center coordinates)
        self.xwake = np.arange(0,self.lx+1e-10,self.dx)
        self.ywake = np.arange(0,self.ly+1e-10,self.dx)
        self.x = np.concatenate((stretch(-self.Lxw,self.Lxw,nx1,rx1)[:-1], self.xwake[:-1], stretch(self.lx,self.Lxe,nx3,rx3)))
        self.y = np.concatenate((stretch(-self.Lys,self.Lys,ny1,ry1)[:-1], self.ywake[:-1], stretch(self.ly,self.Lyn,ny3,ry3)))

        # Blocks
        l(0,'blocks')
        l(0,'(')
        l(1,'hex ( 0 1 2 3 4 5 6 7 )  ( %d %d %d ) simpleGrading ( %.8f %.8f %.8f ) // 0 '%(nx1,ny1,nz1,Rw,Rs,Rz))
        l(1,'hex ( 1 8 9 2 5 10 11 6 )  ( %d %d %d ) simpleGrading ( %.8f %.8f %.8f ) // 0 '%(nx2,ny1,nz1,1.0,Rs,Rz))
        l(1,'hex ( 8 12 13 9 10 14 15 11 )  ( %d %d %d ) simpleGrading ( %.8f %.8f %.8f ) // 0 '%(nx3,ny1,nz1,Re,Rs,Rz))
        l(1,'hex ( 3 2 16 17 7 6 18 19 )  ( %d %d %d ) simpleGrading ( %.8f %.8f %.8f ) // 0 '%(nx1,ny2,nz1,Rw,1.0,Rz))
        l(1,'hex ( 2 9 20 16 6 11 21 18 )  ( %d %d %d ) simpleGrading ( %.8f %.8f %.8f ) // 0 '%(nx2,ny2,nz1,1.0,1.0,Rz))
        l(1,'hex ( 9 13 22 20 11 15 23 21 )  ( %d %d %d ) simpleGrading ( %.8f %.8f %.8f ) // 0 '%(nx3,ny2,nz1,Re,1.0,Rz))
        l(1,'hex ( 17 16 24 25 19 18 26 27 )  ( %d %d %d ) simpleGrading ( %.8f %.8f %.8f ) // 0 '%(nx1,ny3,nz1,Rw,Rn,Rz))
        l(1,'hex ( 16 20 28 24 18 21 29 26 )  ( %d %d %d ) simpleGrading ( %.8f %.8f %.8f ) // 0 '%(nx2,ny3,nz1,1.0,Rn,Rz))
        l(1,'hex ( 20 22 30 28 21 23 31 29 )  ( %d %d %d ) simpleGrading ( %.8f %.8f %.8f ) // 0 '%(nx3,ny3,nz1,Re,Rn,Rz))
        l(0,');')
        l(0,'')

        # Boundaries
        l(0,'boundary')
        l(0,'(')
        l(1,'bottomAndTop')
        l(1,'{')
        l(2,'type empty;')
        l(2,'faces')
        l(2,'(')
        l(3,'(0 1 2 3)')
        l(3,'(4 5 6 7)')
        l(3,'(1 8 9 2)')
        l(3,'(5 10 11 6)')
        l(3,'(8 12 13 9)')
        l(3,'(10 14 15 11)')
        l(3,'(3 2 16 17)')
        l(3,'(7 6 18 19)')
        l(3,'(2 9 20 16)')
        l(3,'(6 11 21 18)')
        l(3,'(9 13 22 20)')
        l(3,'(11 15 23 21)')
        l(3,'(17 16 24 25)')
        l(3,'(19 18 26 27)')
        l(3,'(16 20 28 24)')
        l(3,'(18 21 29 26)')
        l(3,'(20 22 30 28)')
        l(3,'(21 23 31 29)')
        l(2,');')
        l(1,'}')

        l(1,'sides')
        l(1,'{')
        l(2,'type patch;')
        l(2,'faces')
        l(2,'(')
        l(3,'(4 5 1 0)')
        l(3,'(5 10 8 1)')
        l(3,'(10 14 12 8)')
        l(3,'(27 26 24 25)')
        l(3,'(26 29 28 24)')
        l(3,'(29 31 30 28)')
        l(2,');')
        l(1,'}')

        l(1,'inlet')
        l(1,'{')
        l(2,'type patch;')
        l(2,'faces')
        l(2,'(')
        l(3,'(4 0 3 7)')
        l(3,'(7 3 17 19)')
        l(3,'(19 17 25 27)')
        l(2,');')
        l(1,'}')

        l(1,'outlet')
        l(1,'{')
        l(2,'type patch;')
        l(2,'faces')
        l(2,'(')
        l(3,'(14 12 13 15)')
        l(3,'(15 13 22 23)')
        l(3,'(23 22 30 31)')
        l(2,');')
        l(1,'}')

        l(0,');')
        l(0,'')

        # Merge patch pairs
        l(0,'mergePatchPairs')
        l(0,'(')
        l(0,');')
        l(0,'')

        # Write the mesh to a blockMeshDict file for OpenFOAM
        self.write_dict2file(self.code_block,filename,'blockMeshDict')

        # Print some statistics
        totalCells = (nx1+nx2+nx3)*(ny1+ny2+ny3)*nz1
        wakeCells = nx2*ny2*nz1
        hypoCells = (self.Lxw+self.lx+self.Lxe)*(self.Lys+self.ly+self.Lyn)/self.dx**2
        print('Total number of cells: %d'%(totalCells))
        #print('Percent of cells in wake domain: %.1f percent'%(wakeCells/totalCells*100))
        print('Cells if uniform grid had been used: %d'%(hypoCells))

    def plot_grid(self,every=4,save_fig=False):
        # Plot the whole domain
        figh = plt.figure()
        gc = 'k'
        wc = 'b'
        lw = 0.5
        # Block bounds
        plt.plot([-self.Lxw, self.lx+self.Lxe],[-self.Lys, -self.Lys],color=gc)
        plt.plot([-self.Lxw, self.lx+self.Lxe],[self.ly+self.Lyn, self.ly+self.Lyn],color=gc)
        plt.plot([-self.Lxw, -self.Lxw],[-self.Lys, self.ly+self.Lyn],color=gc)
        plt.plot([self.lx+self.Lxe, self.lx+self.Lxe],[-self.Lys, self.ly+self.Lyn],color=gc)
        # Grid lines
        xwhole = self.x[::every]
        ywhole = self.y[::every]
        for i in range(len(xwhole)):
            plt.plot([xwhole[i],xwhole[i]], [-self.Lys, self.ly+self.Lyn], color=gc, linewidth=lw)
        for i in range(len(ywhole)):
            plt.plot([-self.Lxw,self.lx+self.Lxe], [ywhole[i], ywhole[i]], color=gc, linewidth=lw)
        # Plot whole domain
        plt.plot([0, self.lx],[0, 0],color=wc)
        plt.plot([0, self.lx],[self.ly, self.ly],color=wc)
        plt.plot([0, 0],[0, self.ly],color=wc)
        plt.plot([self.lx, self.lx],[0, self.ly],color=wc)
        # Turbines
        def draw_turbine(wt):
            plt.plot([wt['x'], wt['x']],[wt['y']-wt['D']/2,wt['y']+wt['D']/2],color='green')
        for i in range(len(self.wts)):
            draw_turbine(self.wts[i])
        # to be done
        plt.axis('scaled')
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
        plt.title('Whole domain (every %dth grid line)'%every)

        # Plot the wake domain
        figw = plt.figure()
        # Outer bound
        plt.plot([0, self.lx],[0, 0],color=wc)
        plt.plot([0, self.lx],[self.ly, self.ly],color=wc)
        plt.plot([0, 0],[0, self.ly],color=wc)
        plt.plot([self.lx, self.lx],[0, self.ly],color=wc)
        # Grid lines
        xwake = self.xwake[::every]
        ywake = self.ywake[::every]
        for i in range(len(xwake)):
            plt.plot([xwake[i],xwake[i]], [0, self.ly], color=gc, linewidth=lw)
        for i in range(len(ywake)):
            plt.plot([0,self.lx], [ywake[i], ywake[i]], color=gc, linewidth=lw)
        # Turbines
        for i in range(len(self.wts)):
            draw_turbine(self.wts[i])
        plt.axis('scaled')
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
        plt.title('Wake domain (every %dth grid line shown)'%every)

        if save_fig == True:
            pdf = mpl.backends.backend_pdf.PdfPages("grid.pdf")
            pdf.savefig(figh, bbox_inches='tight')
            pdf.savefig(figw, bbox_inches='tight')
            pdf.close()



    def make_topo(self,filename):
        '''
            Create blockMeshDict
        '''

        # String to be made
        self.code_topo = ''
        def l(id,s):
            self.code_topo += lineOut(id,s)

        # Write a single turbine to topoSetDict
        def write_turbine(n,wt):
            l(1,'// actuationDisk%d'%n)
            l(1,'{')
            l(2,'name    actuationDisk%dCellSet;'%n)
            l(2,'type    cellSet;')
            l(2,'action  new;')
            l(2,'source  cylinderToCell;')
            l(2,'sourceInfo')
            l(2,'{')
            l(3,'p1     (%.1f %.1f %.1f);'%(wt['x']-self.dx, wt['y'], self.lz/2))
            l(3,'p2     (%.1f %.1f %.1f);'%(wt['x']+self.dx, wt['y'], self.lz/2))
            l(3,'radius %.1f;'%(wt['D']/2))
            l(2,'}')
            l(1,'}')
            l(1,'{')
            l(2,'name    actuationDisk%d;'%n)
            l(2,'type    cellZoneSet;')
            l(2,'action  new;')
            l(2,'source  setToCellZone;')
            l(2,'sourceInfo')
            l(2,'{')
            l(3,'set actuationDisk%dCellSet;'%n)
            l(2,'}')
            l(1,'}')
            l(0,'')

        # Write the topoSetDict string
        l(0,'actions')
        l(0,'(')
        for i in range(len(self.wts)):
            write_turbine(i+1,self.wts[i])
        l(0,');')

        # Write topoSetDict file
        self.write_dict2file(self.code_topo,filename,'topoSetDict')


    def make_fvoptions(self,filename):
        '''
            Create fvOptions.
            Also makes a text file with simple turbine info.
        '''

        # String to be made
        self.code_fvoptions = ''
        def l(id,s):
            self.code_fvoptions += lineOut(id,s)

        def write_turbine(n,wt):
            l(0,'// Disk %d'%n)
            l(0,'disk%d'%n)
            l(0,'{')
            l(1,'type            actuatorDiskFoam;')
            l(1,'variant         calaf;')
            l(1,'selectionMode   cellSet;')
            l(1,'cellSet         actuationDisk%d;'%n)
            l(1,'diskDir         (1 0 0);')
            l(1,'Ct              %.4f;'%wt['CT'])
            l(1,'diskArea        %.1f;'%wt['D'])
            l(0,'}')
            l(0,'')

        # Write the turbines to fvOptions.
        for i in range(len(self.wts)):
            write_turbine(i+1,self.wts[i])

        # Write turbine info to simple text file.
        with open('turbines.dat', "w") as file:
            file.write('Turbine number, x, y, z, D')
            for i in range(len(self.wts)):
                turbi = self.wts[i]
                file.write('\n%d %.4f %.4f %.4f %.4f'%(i+1, turbi['x'], turbi['y'], turbi['zh'], turbi['D']))

        # Write source term to avoid flow development
        if self.kesource == True:
            l(0,'// K- and eps-source terms to avoid flow development')
            l(0,'ScalarSemiImplicitSource1')
            l(0,'{')
            l(1,'type                scalarSemiImplicitSource;')
            l(1,'selectionMode       all;')
            l(1,'volumeMode          specific;')
            l(1,'// Specification of sources in OpenFOAM-2206 and newer)')
            l(1,'sources')
            l(1,'{')
            l(2,'// Specified as ( explicit(Su), implicit(Sp) ):')
            l(2,'k           (%.10f 0);'%(self.Sk))
            l(2,'epsilon     (%.10e  0);'%(self.Seps))
            l(1,'}')
            l(0,'}')

        # Write fvOptions file
        self.write_dict2file(self.code_fvoptions,filename,'fvOptions')

    def make_inflow(self,filename):
        '''
            Create ASLparameters file.
        '''

        # String to be made
        self.code_inflow = ''
        def l(id,s):
            self.code_inflow += lineOut(id,s)

        # Header info
        l(0,'// Paramters for the neutral atmospheric surface layer (ASL)')
        l(0,'// aka. the log-law inflow.')
        l(0,'Uref                 %.4f;'%(self.Uinf))
        l(0,'flowVelocity         ($Uref 0 0);')
        l(0,'pressure             0;')
        l(0,'turbulentKE          %.10f;'%(self.Kinf))
        l(0,'turbulentEpsilon     %.10f;'%(self.epsinf))
        l(0,'')
        l(0,'// These parameters are not used by the code, but were used to calculate the above.')
        l(0,'//Cmu                %.10f'%(self.Cmu) + ' [-]')
        l(0,'//kappa              %.10f'%(self.kap) + ' [-]')
        l(0,'//z0                 %.10f'%(self.z0) + ' [m]')
        l(0,'//utau               %.10f'%(self.utau) + ' [m/s]')
        l(0,'//zref               %.1f'%(self.zref) + ' [m]')
        l(0,'//Iref               %.2f'%(self.Iinf*100) + ' [%]')
        l(0,'')

        # Write ASLparameters file
        self.write_dict2file(self.code_inflow,filename,ff=False)

    def make_turbulence(self,filename):
        '''
            Create turbulenceProperties file.
        '''

        # String to be made
        self.code_turbulence = ''
        def l(id,s):
            self.code_turbulence += lineOut(id,s)

        l(0,'simulationType RAS;')
        l(0,'RAS')
        l(0,'{')
        l(1,'RASModel            kEpsilon;')
        l(1,'turbulence          on;')
        l(1,'printCoeffs         on;')
        l(1,'kEpsilonCoeffs')
        l(1,'{')
        l(2,'Cmu         %.10f;'%self.Cmu)
        l(2,'C1          %.10f;'%self.Ce1)
        l(2,'C2          %.10f;'%self.Ce2)
        l(2,'sigmaEps    %.10f;'%self.sige)
        l(1,'}')
        l(0,'}')
        l(0,'')

        # Write turbulenceProperties file
        self.write_dict2file(self.code_turbulence,filename,ob='turbulenceProperties')



if __name__ == "__main__":
    print('Examples are available in farm2d/examples')

