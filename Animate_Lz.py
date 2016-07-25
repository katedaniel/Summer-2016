import numpy as np
import matplotlib.pyplot as plt
import Orbit_Code as OC
import matplotlib.animation as animation
from matplotlib.colors import LogNorm
from timeit import default_timer
reload(OC)
  
print("Load stuff for plotting...")
start = default_timer()

dpi = 300

theta = "30"
Lz = np.load('/Users/LBarbano/Desktop/QP_Data/Lz_data_(theta='+theta+').npy')
Lz = Lz[:,0:1000]
Lz0 = Lz[:,0][:, np.newaxis]
del_Lz = Lz-Lz0
Lz0 = Lz[:,0]
duration = default_timer() - start
print "time: %s s" % str(duration) 

RC = 1 #Retarding constant for animation

#This function animates a surface density histogram of Lz0 versus delta Lz
def animate_Lz_hist():
    
    plt.close('all')    
    fig = plt.figure() #create figure
    ax = fig.add_subplot(111,aspect = 'equal')  #create and add nx4 subplots
    
    ax.set_title("Change in Angular Momentum ($theta$ = 30) Time(Gy): %f" % 0.,size = 20)
    ax.axvline(1760., color='black',lw = 1) #corotation
    ax.axvline(2382.25,color='black', ls='dashed',lw =2) #outer lindblad
    ax.axvline(1137.746,color='black', ls='dashed',lw =2)#inner lindblad
    ax.axvline(2071.13,color='black', ls='dotted', lw=3) #outer ultraharmonic
    ax.axvline(1448.87,color='black', ls='dotted', lw=3) #inner ultraharmonic 
    fig.text(0.5, 0.02, r'Initial $L (1000 kpc\frac{km}{s})$', ha='center', size=20)
    fig.text(0.02, 0.5, r'$\Delta L (1000 kpc\frac{km}{s})$', va='center', rotation='vertical', size=20)
    fig.text(0.85, 0.87, 'Number of\nStars (log)') 
    
    my_cmap = plt.cm.get_cmap('spectral')
    H,xedges,yedges = np.histogram2d(Lz[:,0], del_Lz[:,500], range=[[400,3500], [-1000,1000]],bins=(100, 80)) 
    im = ax.imshow(H.T,origin='low',extent=[700,3000,-1000,1000],interpolation='nearest',aspect='auto',cmap=my_cmap,norm=LogNorm(),vmin=1.)
    plt.colorbar(im)
    anim = animation.FuncAnimation(fig, update_hist, frames= len(Lz.transpose())/RC, interval= 1.0, 
    fargs=(fig,ax,im))
    plt.show()
    
    writer = animation.writers['ffmpeg'](fps=30)
    start = default_timer()
    print("Creating hist movie...")
    anim.save('/Users/LBarbano/Desktop/QP_Data/Animate_Lz_hist(theta='+theta+').mp4',writer=writer,dpi=dpi)
    duration = default_timer() - start
    print "time: %s s" % str(duration) 
    return anim
       
#This function animates a scatter plot of Lz0 versus delta Lz and colors the data points
#according to the surface density histogram of Lz0 versus delta Lz         
def animate_Lz_scat():
    x = Lz0
    y = del_Lz[:,500]
    
    #Figure stuff
    plt.close('all')    
    fig = plt.figure() #create figure
    ax = fig.add_subplot(111,aspect = 'equal')  #create and add nx4 subplots
    ax.set_title("Change in Angular Momentum ($theta$ = 30) Time(Gy): %f" % 0.)
    ax.axvline(1760., color='black',lw = 1) #corotation
    ax.axvline(2382.25,color='black', ls='dashed',lw =2) #outer lindblad
    ax.axvline(1137.746,color='black', ls='dashed',lw =2)#inner lindblad
    ax.axvline(2071.13,color='black', ls='dotted', lw=3) #outer ultraharmonic
    ax.axvline(1448.87,color='black', ls='dotted', lw=3) #inner ultraharmonic 
    fig.text(0.5, 0.02, r'Initial $L (1000 kpc\frac{km}{s})$', ha='center', size=20)
    fig.text(0.02, 0.5, r'$\Delta L (1000 kpc\frac{km}{s})$', va='center', rotation='vertical', size=20)
    fig.text(0.85, 0.87, 'Number of\nStars (log)')     
    
    #Histogram stuff
    my_cmap = plt.cm.get_cmap('spectral')
    H,xedges,yedges = np.histogram2d(x,y, range=[[700,3000], [-1000,1000]],bins=(80, 60)) 
    xidx = np.clip(np.digitize(x, xedges), 0, H.shape[0]-1)
    yidx = np.clip(np.digitize(y, yedges), 0, H.shape[1]-1)
    col = H[xidx, yidx] +1
    scat = ax.scatter(Lz[:,0],del_Lz[:,0],c = col, s = 20,cmap=my_cmap,
    vmax = col.max(),norm=LogNorm(),edgecolor = 'none') #Draw initial location of star 
    plt.colorbar(scat)
    plt.axis([700,3000,-1000,1000])
    
    #Animation stuff
    anim = animation.FuncAnimation(fig, update_scat, frames= len(Lz.transpose())/(2*RC), interval= 1, 
    fargs=(fig,ax,scat))
    writer = animation.writers['ffmpeg'](fps=30)
    start = default_timer()
    print("Creating scat movie...")
    anim.save('/Users/LBarbano/Desktop/QP_Data/Animate_Lz_scat(theta='+theta+').mp4',writer=writer,dpi=dpi)
    duration = default_timer() - start
    print "time: %s s" % str(duration) 
    return anim
    
    
def update_scat(i,fig,ax,scat):
    x = Lz0
    y = del_Lz[:,i*RC*2]
    data = np.hstack((x[:,np.newaxis], y[:, np.newaxis]))
    scat.set_offsets(data)    
    H,xedges,yedges = np.histogram2d(x,y, range=[[700,3500], [-1000,1000]],bins=(100,80)) 
    xidx = np.clip(np.digitize(x, xedges), 0, H.shape[0]-1)
    yidx = np.clip(np.digitize(y, yedges), 0, H.shape[1]-1)
    col = H[xidx, yidx] +1
    scat.set_array(col)
    ax.set_title("Change in Angular Momentum ("+ r'$\theta='+theta+'$'+") \n Time(Gy): %f" % float(2*i*RC/1000.),size = 25)
    return scat
    
def update_hist(i,fig,ax,im):
    H,x,y = np.histogram2d(Lz[:,0], del_Lz[:,i*RC], range=[[400,3500], [-1000,1000]],bins=(100, 80)) 
    im.set_array(H.T) 
    ax.set_title("Change in Angular Momentum ("+ r'$\theta='+theta+'$'+") \n Time(Gy): %f" % float(i*RC/1000.),size = 20)
    return im,

    