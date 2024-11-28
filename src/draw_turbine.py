import matplotlib.pyplot as plt
import numpy as np


def draw_turbine(D = 10):
    c = 'k'
    hub_r = D/40
    hub_l = D/10
    root_l = D/10
    root_r = D/20
    def generate_semicircle(center_x, center_y, radius, stepsize=0.1):
        """
        generates coordinates for a semicircle, centered at center_x, center_y
        """        
    
        x = np.arange(center_x, center_x+radius+stepsize, stepsize)
        y = np.sqrt(radius**2 - x**2)
    
        # since each x value has two corresponding y-values, duplicate x-axis.
        # [::-1] is required to have the correct order of elements for plt.plot. 
        x = np.concatenate([x,x[::-1]])
    
        # concatenate y and flipped y. 
        y = np.concatenate([y,-y[::-1]])
    
        return -x, y + center_y
    
    def generate_wing(r,beta):
        
        x = beta/r
        y = r
        
        return x, y
    
    
    
    beta = root_l*(hub_r+root_r)
    tip_l = beta/(D/2)
    
    rtop = np.linspace(hub_r + root_r,D/2,100)
    
    xtmp,ytmp = generate_semicircle(0, 0, hub_r, stepsize=D/200)
    x1,y1 = generate_wing(rtop,beta)
    
    
    # Hub
    plt.plot(xtmp,ytmp,color=c)
    plt.plot([0,hub_l],[-hub_r,-hub_r],color=c)
    plt.plot([0,hub_l],[hub_r,hub_r],color=c)
    plt.plot([hub_l,hub_l],[-hub_r,hub_r],color=c)
    
    # Upper wing
    plt.plot([0,0],[hub_r,D/2],color=c)
    plt.plot([0,tip_l],[D/2,D/2],color=c)
    plt.plot(x1,y1,color=c)
    plt.plot([0,root_l],[hub_r,hub_r+root_r],color=c)
    
    # Lower wing
    plt.plot([0,0],[-D/2,-hub_r],color=c)
    plt.plot([0,tip_l],[-D/2,-D/2],color=c)
    plt.plot(x1,-y1,color=c)
    plt.plot([0,root_l],[-hub_r,-hub_r-root_r],color=c)



if __name__ == "__main__":
    D = 10
    fig, ax = plt.subplots(figsize=(6, 6))
    draw_turbine(D)
    plt.axis('scaled')
    plt.xlim([-D,D])
    plt.ylim([-D,D])