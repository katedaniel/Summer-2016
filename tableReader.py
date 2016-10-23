import os
import numpy as np
import matplotlib.pyplot as plt
import Orbit_Code
from matplotlib.colors import LogNorm
reload(Orbit_Code)

filepath = "C:/Users/Noah/Documents/GitHub/Trapped_Orbital_Integrator/" 

#Get names of text files        
files1 = [i for i in os.listdir(filepath) if os.path.isfile(os.path.join(filepath,i)) and 'table1' in i]
files2 = [i for i in os.listdir(filepath) if os.path.isfile(os.path.join(filepath,i)) and 'table2' in i]
files_ROE = [i for i in os.listdir(filepath) if os.path.isfile(os.path.join(filepath,i)) and 'table_ROE' in i]

#Import all data
tableInfo1 = np.array([np.loadtxt(filepath+files1[i],delimiter=" ",dtype= str)[:,1:16].astype(float) for i in range(len(files1))]) 
tableInfo2 = np.array([np.loadtxt(filepath+files2[i],delimiter=" ",dtype= str).astype(float) for i in range(len(files2))]) 
tableInfo_ROE = np.array([np.loadtxt(filepath+files_ROE[i],delimiter=" ",dtype= str).astype(float) for i in range(len(files_ROE))]) 

###This function plots fraction of trapped stars over time for each theta
def trap_plot():
    
    t = tableInfo2[0][:,0] 
    trap_frac = np.array([tableInfo2[i][:,1] for i in range(len(tableInfo2))])
    lam_spec = np.array([tableInfo1[i][:,9] for i in range(len(tableInfo1))])
    start_trapped = [(lam_spec[i] == 0).sum()+(lam_spec[i] == 1).sum()+(lam_spec[i] == 2).sum() for i in range(len(tableInfo2))]
    
    plt.close('all')
    colors = ['green','blue','purple','black']
    labels = ['Theta = 15 (%s ','Theta = 20 (%s ','Theta = 25 (%s ','Theta = 30 (%s ']
    [plt.plot(t, trap_frac[i]/trap_frac[i][0], label=labels[i] % (start_trapped[i]) +'trapped stars at t=0)', color=colors[i])for i in range(len(tableInfo2))]
    
    plt.xlabel('Time (years)',size=18)
    plt.ylabel('Trapped stars (normalized)',size=18)
    plt.title('Fraction of Trapped Stars per Theta', size=22)
    plt.legend()
    plt.show()

###This function is basically the previous trap plot, but seperated by resonance and spiral arm interactions   
def trap_plot2():
    t = tableInfo2[0][:,0] 
    trap_frac = np.array([tableInfo2[i][:,1] for i in range(len(tableInfo2))])
    trap_frac1 = np.array([tableInfo2[i][:,4] for i in range(len(tableInfo2))])
    trap_frac2 = np.array([tableInfo2[i][:,5] for i in range(len(tableInfo2))])

    plt.close('all')
    fig = plt.figure() #create figure
    fig.subplots_adjust(wspace=0.03,hspace=0.03) #some plotting stuff
    ax = [fig.add_subplot(2,2,i+1,aspect = 'auto') for i in range(4)] #create and add 2x2 subplots
    
    for i in range(4):
        ax[i].plot(t,trap_frac[i]/trap_frac[i][0],color='black')
        ax[i].plot(t,trap_frac1[i]/trap_frac[i][0],color='blue')
        ax[i].plot(t,trap_frac2[i]/trap_frac[i][0],color='red')
    
    ax[0].set_xticklabels([])
    ax[1].set_xticklabels([])
    ax[1].set_yticklabels([])
    ax[3].set_yticklabels([])
    ax[1].legend(['Total','No Interaction','Resonance Interaction'])
    plt.suptitle('Fraction of Trapped Stars', size=22)
    fig.text(0.5, 0.02, 'Time (years)', ha='center', size=18)
    fig.text(0.02, 0.5, 'Number of Trapped Stars (normalized)', va='center', rotation='vertical', size=18)
    
    plt.show()
    
###This function plots rms of change in angmom over time for each theta
def angmom_plot():
      
    t = tableInfo2[0][:,0] 
    Lz_rms = np.array([tableInfo2[i][:,2] for i in range(len(tableInfo2))])
    Lz_rms_spec = np.array([tableInfo2[i][:,3] for i in range(len(tableInfo2))])
    lam_spec = np.array([tableInfo1[i][:,9] for i in range(len(tableInfo1))])
    start_trapped = [(lam_spec[i] == 0).sum()+(lam_spec[i] == 1).sum()+(lam_spec[i] == 2).sum() for i in range(len(tableInfo1))]
    
    plt.close('all')
    colors = ['green','blue','purple','black']
    labels = ['Theta = 15','Theta = 20','Theta = 25','Theta = 30']
    [plt.plot(t, Lz_rms[i]/start_trapped[i], label=labels[i], color=colors[i])for i in range(len(tableInfo2))]
    [plt.plot(t, Lz_rms_spec[i]/start_trapped[i], color=colors[i], ls='dashed')for i in range(len(tableInfo2))]
    
    plt.xlabel('Time (years)',size=18)
    plt.ylabel(r'rms of $\Delta$L ($kpc\frac{km}{s}$) (normalized)', size=18)
    plt.title('Change in Angular Momenta per Theta', size=22)
    plt.legend(loc='upper left')
    plt.show()
    
###This function returns rms of angmom for each theta at different time intervals
def Lz_rms():
    num = [11,12,13,14]
    rms = np.array([np.sqrt((tableInfo1[i][:,[10,num[j]]]**2).mean(axis=1)) for i in range(len(tableInfo1)) for j in range(len(num))])
    return np.split(rms,4)

###This function returns the fraction of initially trapped orbits that end trapped for each theta
def trap_frac():
    lam_spec = np.array([tableInfo1[i][:,9] for i in range(len(tableInfo1))])
    start_trapped = np.array([float((lam_spec[i] == 0).sum()+(lam_spec[i] == 1).sum()+(lam_spec[i] == 2).sum()) for i in range(len(tableInfo1))])
    end_trapped = np.array([float((lam_spec[i] == 0).sum()+(lam_spec[i] == 1).sum()) for i in range(len(tableInfo1))])
    #calculate fraction of initially trapped orbits that end trapped
    trap_frac = end_trapped/start_trapped
    return trap_frac

def plot_Lz():
    #Extract Lz data and calculate delta Lz
    Lz = np.array([tableInfo1[i][:,10:15] for i in range(len(tableInfo1))])
    del_Lz = np.array([np.subtract(Lz[i][:,j+1],Lz[i][:,0]) for i in range(len(tableInfo1)) for j in range (len(Lz[i].transpose())-1)])

    #It's plotting time baby ohhhhh yeah leggoooo
    plt.close('all')    
    
    dimx = len(files1) #n rows of subplots (can be any number depending of how many thetas we want, in this case 4)
    fig = plt.figure() #create figure
    fig.subplots_adjust(wspace=0,hspace=0) #some plotting stuff
    ax = [fig.add_subplot(dimx,4,i+1,aspect = 'auto') for i in range(dimx*4)] #create and add nx4 subplots
    
    #Generate histogram stuff
    histograms = np.array([np.histogram2d(Lz[i/4][:,0], del_Lz[i], range=[[500,3500], [-500,500]],bins=(80, 60)) for i in range(len(del_Lz)) ])
    my_cmap = plt.cm.get_cmap('spectral')
    my_cmap.set_under('w')
    im = [ax[i].imshow(histograms[i,0].T,origin='low',extent=[500.0, 3500.0, -500.0, 500.0],interpolation='nearest',aspect='auto', cmap=my_cmap,norm=LogNorm(),vmin=1.) for i in range(len(ax))]
    
    #Hardcode some plot adjustments
    #the array lists in the below for loops contain the indeces of the subplots to be adjusted
    [ax[i-1].set_yticklabels([])  for i in [2,3,4,6,7,8,10,11,12,14,15,16]]
    [ax[i-1].set_xticklabels([])  for i in [1,5,9,4,8,12]]
    [ax[i-1].set_xticklabels(['',1,'',2,'',3])  for i in [13,14,15,16]]
    [ax[i-1].set_yticklabels(['',-0.4,-0.2,0,0.2,0.4])  for i in [1,5,9,13]]
    
    theta_labels = [r'$\theta=15$',r'$\theta=20$',r'$\theta=25$',r'$\theta=30$']
    time_labels = ['t = 0.5','t = 1.0','t = 1.5','t = 2.0']
    [ax[val-1].set_ylabel(theta_labels[i])  for i, val in enumerate([1,5,9,13])]
    [ax[val-1].set_xlabel(time_labels[i])  for i, val in enumerate([13,14,15,16])]
        
    #Plot corotation and lindblad resonances
    
    for i in range(len(ax)):
        ax[i].axvline(1760., color='black') #corotation
        ax[i].axvline(2382.25,color='black', ls='dashed') #outer lindblad
        ax[i].axvline(1137.746,color='black', ls='dashed')#inner lindblad
        ax[i].axvline(2071.13,color='black', ls='dotted', lw=2) #outer ultraharmonic
        ax[i].axvline(1448.87,color='black', ls='dotted', lw=2) #inner ultraharmonic        

    #adding overall axis labels, title, and colorbar label
    fig.text(0.5, 0.02, r'Initial L $(1000 kpc\frac{km}{s})$', ha='center', size=18)
    fig.text(0.02, 0.5, r'$\Delta L (1000 kpc\frac{km}{s})$', va='center', rotation='vertical', size=18)
    fig.text(0.85, 0.87, 'Number of\nStars (log)')
    fig.text(0.905,0.75,'40')
    fig.text(0.905,0.63,'20')
    fig.text(0.905,0.5,'10')
    fig.text(0.905,0.41,'5')
    fig.text(0.905,0.26,'2')
    fig.text(0.905,0.14,'1')
    plt.suptitle('Change in Angular Momentum', size=22)
    
    #slight adjustments and put it all together
    fig.subplots_adjust(right=0.8, bottom=0.17, left=0.17)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im[0],cax=cbar_ax, ticks=[0.9,2,5,9,20,40])
    plt.show()
    
def ROE_hist():
    #Extract ROE data ad del_ROE
    ROE_tot15 = np.array(tableInfo_ROE[0])
    ROE_tot20 = np.array(tableInfo_ROE[1])
    ROE_tot25 = np.array(tableInfo_ROE[2])
    ROE_tot30 = np.array(tableInfo_ROE[3])
    ROE15_1 = []
    ROE20_1 = []
    ROE25_1 = []
    ROE30_1 = []
    ROE15_2 = []
    ROE20_2 = []
    ROE25_2 = []
    ROE30_2 = []
    for ROE in ROE_tot15:
        if ROE[0] == 1.:
            ROE15_1.append(ROE[1:])
        if ROE[0] == 2.:
            ROE15_2.append(ROE[1:])
    for ROE in ROE_tot20:
        if ROE[0] == 1.:
            ROE20_1.append(ROE[1:])
        if ROE[0] == 2.:
            ROE20_2.append(ROE[1:])
    for ROE in ROE_tot25:
        if ROE[0] == 1.:
            ROE25_1.append(ROE[1:])
        if ROE[0] == 2.:
            ROE25_2.append(ROE[1:])
    for ROE in ROE_tot30:
        if ROE[0] == 1.:
            ROE30_1.append(ROE[1:])
        if ROE[0] == 2.:
            ROE30_2.append(ROE[1:])
    
    del_ROE = np.array([np.subtract(ROE[i][:,j+1],ROE[i][:,0]) for i in range(len(tableInfo_ROE)) for j in range (len(ROE[i].transpose())-1)])
    
    #It's plotting time baby ohhhhh yeah leggoooo
    plt.close('all')    
    
    fig = plt.figure() #create figure
    fig.subplots_adjust(wspace=0,hspace=0) #some plotting stuff
    ax = [fig.add_subplot(4,2,i+1,aspect = 'auto') for i in range(8)] #create and add nx2 subplots
    
    #Generate histogram stuff
    histogram15_1 = np.histogram2d(ROE15_1, del_Lz[i], range=[[500,3500], [-500,500]],bins=(80, 60))
    
    histograms = np.array([np.histogram2d(ROE[:,0], del_Lz[i], range=[[500,3500], [-500,500]],bins=(80, 60)) for i in range(len(del_Lz))])
    my_cmap = plt.cm.get_cmap('spectral')
    my_cmap.set_under('w')
    im = [ax[i].imshow(histograms[i,0].T,origin='low',extent=[500.0, 3500.0, -500.0, 500.0],interpolation='nearest',aspect='auto', cmap=my_cmap,norm=LogNorm(),vmin=1.) for i in range(len(ax))]
    