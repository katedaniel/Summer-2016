import numpy as np
import matplotlib.pyplot as plt

#function that finds root mean square
def rmse(values):
    return np.sqrt((values**2).mean())

dataFilePath_25 = "C:/Users/Noah/Documents/GitHub/Trapped_Orbital_Integrator/table_25.txt"
tableInfo_25 = np.loadtxt(dataFilePath_25,delimiter=" ",dtype= str)
filepaths_25 = tableInfo_25[:,0]
data_25 = tableInfo_25[:,1:16].astype(float)
dataFilePath_15 = "C:/Users/Noah/Documents/GitHub/Trapped_Orbital_Integrator/table_15.txt"
tableInfo_15 = np.loadtxt(dataFilePath_15,delimiter=" ",dtype= str)
filepaths_15 = tableInfo_15[:,0]
data_15 = tableInfo_15[:,1:16].astype(float)
dataFilePath_20 = "C:/Users/Noah/Documents/GitHub/Trapped_Orbital_Integrator/table_20.txt"
tableInfo_20 = np.loadtxt(dataFilePath_20,delimiter=" ",dtype= str)
filepaths_20 = tableInfo_20[:,0]
data_20 = tableInfo_20[:,1:16].astype(float)
dataFilePath_30 = "C:/Users/Noah/Documents/GitHub/Trapped_Orbital_Integrator/table_30.txt"
tableInfo_30 = np.loadtxt(dataFilePath_30,delimiter=" ",dtype= str)
filepaths_30 = tableInfo_30[:,0]
data_30 = tableInfo_30[:,1:16].astype(float)
'''
rms1 = np.sqrt((data[:,[10,11]]**2).mean(axis=1)) #rms 0 to 0.5
rms2 = np.sqrt((data[:,[10,12]]**2).mean(axis=1)) #rms 0 to 1
rms3 = np.sqrt((data[:,[10,13]]**2).mean(axis=1)) #rms 0 to 1.5
rms4 = np.sqrt((data[:,[10,14]]**2).mean(axis=1)) #rms 0 to 2.0

lam_spec = data[:,9]
start_trapped = float(((lam_spec == 0) or (lam_spec == 1) or (lam_spec == 2)).sum())
end_trapped = float(((lam_spec == 0) or (lam_spec == 1)).sum())
trap_frac = end_trapped/start_trapped
'''
def Lz_plot():
    #pulling angmom data
    Lz_15 = data_15[:,[10,11,12,13,14]]
    Lz_20 = data_20[:,[10,11,12,13,14]]
    Lz_25 = data_25[:,[10,11,12,13,14]]
    Lz_30 = data_30[:,[10,11,12,13,14]]
    #data analysis for THETA 15
    Lz_15_0 = Lz_15[:,0]#initial angmom
    Lz_15_1 = Lz_15[:,1]#angmom at t=0.5
    Lz_15_2 = Lz_15[:,2]#angmom at t=1.0
    Lz_15_3 = Lz_15[:,3]#angmom at t=1.5
    Lz_15_4 = Lz_15[:,4]#angmom at t=2.0
    del_15_1 = np.subtract(Lz_15_1,Lz_15_0)#change in angmom after t=0.5
    del_15_2 = np.subtract(Lz_15_2,Lz_15_0)#change in angmom after t=1.0
    del_15_3 = np.subtract(Lz_15_3,Lz_15_0)#change in angmom after t=1.5
    del_15_4 = np.subtract(Lz_15_4,Lz_15_0)#change in angmom after t=2.0
    #data analysis for THETA 20
    Lz_20_0 = Lz_20[:,0]
    Lz_20_1 = Lz_20[:,1]
    Lz_20_2 = Lz_20[:,2]
    Lz_20_3 = Lz_20[:,3]
    Lz_20_4 = Lz_20[:,4]
    del_20_1 = np.subtract(Lz_20_1,Lz_20_0)
    del_20_2 = np.subtract(Lz_20_2,Lz_20_0)
    del_20_3 = np.subtract(Lz_20_3,Lz_20_0)
    del_20_4 = np.subtract(Lz_20_4,Lz_20_0)
    #data analysis for THETA 15
    Lz_25_0 = Lz_25[:,0]
    Lz_25_1 = Lz_25[:,1]
    Lz_25_2 = Lz_25[:,2]
    Lz_25_3 = Lz_25[:,3]
    Lz_25_4 = Lz_25[:,4]
    del_25_1 = np.subtract(Lz_25_1,Lz_25_0)
    del_25_2 = np.subtract(Lz_25_2,Lz_25_0)
    del_25_3 = np.subtract(Lz_25_3,Lz_25_0)
    del_25_4 = np.subtract(Lz_25_4,Lz_25_0)
    #data analysis for THETA 15
    Lz_30_0 = Lz_30[:,0]
    Lz_30_1 = Lz_30[:,1]
    Lz_30_2 = Lz_30[:,2]
    Lz_30_3 = Lz_30[:,3]
    Lz_30_4 = Lz_30[:,4]
    del_30_1 = np.subtract(Lz_30_1,Lz_30_0)
    del_30_2 = np.subtract(Lz_30_2,Lz_30_0)
    del_30_3 = np.subtract(Lz_30_3,Lz_30_0)
    del_30_4 = np.subtract(Lz_30_4,Lz_30_0)
    
    plt.close('all')
    
    #MAKING THE SUBPLOTS
    
    fig = plt.figure()
    fig.subplots_adjust(wspace=0,hspace=0)
    color = 'OrRd'
    
    #THETA OF 15
    #t=0.5
    ax1 = fig.add_subplot(4,4,1,aspect='auto')
    H1, xedges1, yedges1 = np.histogram2d(Lz_15_0, del_15_1, range=[[500,3500], [-500,500]], bins=(80, 60))
    myextent1  =[xedges1[0],xedges1[-1],yedges1[0],yedges1[-1]]
    im = plt.imshow(H1.T,origin='low',extent=myextent1,interpolation='nearest',aspect='auto', cmap=color)
    #t=1.0
    ax2 = fig.add_subplot(4,4,2,aspect='auto')
    H2, xedges2, yedges2 = np.histogram2d(Lz_15_0, del_15_2, range=[[500,3500], [-500,500]], bins=(80, 60))
    myextent2  =[xedges2[0],xedges2[-1],yedges2[0],yedges2[-1]]
    plt.imshow(H2.T,origin='low',extent=myextent2,interpolation='nearest',aspect='auto', cmap=color)
    #t=1.5
    ax3 = fig.add_subplot(4,4,3,aspect='auto')
    H3, xedges3, yedges3 = np.histogram2d(Lz_15_0, del_15_3, range=[[500,3500], [-500,500]], bins=(80, 60))
    myextent3  =[xedges3[0],xedges3[-1],yedges3[0],yedges3[-1]]
    plt.imshow(H3.T,origin='low',extent=myextent3,interpolation='nearest',aspect='auto', cmap=color)
    #t=2.0
    ax4 = fig.add_subplot(4,4,4,aspect='auto')
    H4, xedges4, yedges4 = np.histogram2d(Lz_15_0, del_15_4, range=[[500,3500], [-500,500]], bins=(80, 60))
    myextent4  =[xedges4[0],xedges4[-1],yedges4[0],yedges4[-1]]
    plt.imshow(H4.T,origin='low',extent=myextent4,interpolation='nearest',aspect='auto', cmap=color)
    
    #THETA OF 20
    #t=0.5
    ax5 = fig.add_subplot(4,4,5,aspect='auto')
    H1, xedges1, yedges1 = np.histogram2d(Lz_20_0, del_20_1, range=[[500,3500], [-500,500]], bins=(80, 60))
    myextent1  =[xedges1[0],xedges1[-1],yedges1[0],yedges1[-1]]
    plt.imshow(H1.T,origin='low',extent=myextent1,interpolation='nearest',aspect='auto', cmap=color)
    #t=1.0
    ax6 = fig.add_subplot(4,4,6,aspect='auto')
    H2, xedges2, yedges2 = np.histogram2d(Lz_20_0, del_20_2, range=[[500,3500], [-500,500]], bins=(80, 60))
    myextent2  =[xedges2[0],xedges2[-1],yedges2[0],yedges2[-1]]
    plt.imshow(H2.T,origin='low',extent=myextent2,interpolation='nearest',aspect='auto', cmap=color)
    #t=1.5
    ax7 = fig.add_subplot(4,4,7,aspect='auto')
    H3, xedges3, yedges3 = np.histogram2d(Lz_20_0, del_20_3, range=[[500,3500], [-500,500]], bins=(80, 60))
    myextent3  =[xedges3[0],xedges3[-1],yedges3[0],yedges3[-1]]
    plt.imshow(H3.T,origin='low',extent=myextent3,interpolation='nearest',aspect='auto', cmap=color)
    #t=2.0
    ax8 = fig.add_subplot(4,4,8,aspect='auto')
    H4, xedges4, yedges4 = np.histogram2d(Lz_20_0, del_20_4, range=[[500,3500], [-500,500]], bins=(80, 60))
    myextent4  =[xedges4[0],xedges4[-1],yedges4[0],yedges4[-1]]
    plt.imshow(H4.T,origin='low',extent=myextent4,interpolation='nearest',aspect='auto', cmap=color)
    
    #THETA OF 25
    #t=0.5
    ax9 = fig.add_subplot(4,4,9,aspect='auto')
    H1, xedges1, yedges1 = np.histogram2d(Lz_25_0, del_25_1, range=[[500,3500], [-500,500]], bins=(80, 60))
    myextent1  =[xedges1[0],xedges1[-1],yedges1[0],yedges1[-1]]
    plt.imshow(H1.T,origin='low',extent=myextent1,interpolation='nearest',aspect='auto', cmap=color)
    #t=1.0
    ax10 = fig.add_subplot(4,4,10,aspect='auto')
    H2, xedges2, yedges2 = np.histogram2d(Lz_25_0, del_25_2, range=[[500,3500], [-500,500]], bins=(80, 60))
    myextent2  =[xedges2[0],xedges2[-1],yedges2[0],yedges2[-1]]
    plt.imshow(H2.T,origin='low',extent=myextent2,interpolation='nearest',aspect='auto', cmap=color)
    #t=1.5
    ax11 = fig.add_subplot(4,4,11,aspect='auto')
    H3, xedges3, yedges3 = np.histogram2d(Lz_25_0, del_25_3, range=[[500,3500], [-500,500]], bins=(80, 60))
    myextent3  =[xedges3[0],xedges3[-1],yedges3[0],yedges3[-1]]
    plt.imshow(H3.T,origin='low',extent=myextent3,interpolation='nearest',aspect='auto', cmap=color)
    #t=2.0
    ax12 = fig.add_subplot(4,4,12,aspect='auto')
    H4, xedges4, yedges4 = np.histogram2d(Lz_25_0, del_25_4, range=[[500,3500], [-500,500]], bins=(80, 60))
    myextent4  =[xedges4[0],xedges4[-1],yedges4[0],yedges4[-1]]
    plt.imshow(H4.T,origin='low',extent=myextent4,interpolation='nearest',aspect='auto', cmap=color)
    
    #THETA OF 30
    #t=0.5
    ax13 = fig.add_subplot(4,4,13,aspect='auto')
    H1, xedges1, yedges1 = np.histogram2d(Lz_30_0, del_30_1, range=[[500,3500], [-500,500]], bins=(80, 60))
    myextent1  =[xedges1[0],xedges1[-1],yedges1[0],yedges1[-1]]
    plt.imshow(H1.T,origin='low',extent=myextent1,interpolation='nearest',aspect='auto', cmap=color)
    #t=1.0
    ax14 = fig.add_subplot(4,4,14,aspect='auto')
    H2, xedges2, yedges2 = np.histogram2d(Lz_30_0, del_30_2, range=[[500,3500], [-500,500]], bins=(80, 60))
    myextent2  =[xedges2[0],xedges2[-1],yedges2[0],yedges2[-1]]
    plt.imshow(H2.T,origin='low',extent=myextent2,interpolation='nearest',aspect='auto', cmap=color)
    #t=1.5
    ax15 = fig.add_subplot(4,4,15,aspect='auto')
    H3, xedges3, yedges3 = np.histogram2d(Lz_30_0, del_30_3, range=[[500,3500], [-500,500]], bins=(80, 60))
    myextent3  =[xedges3[0],xedges3[-1],yedges3[0],yedges3[-1]]
    plt.imshow(H3.T,origin='low',extent=myextent3,interpolation='nearest',aspect='auto', cmap=color)
    #t=2.0
    ax16 = fig.add_subplot(4,4,16,aspect='auto')
    H4, xedges4, yedges4 = np.histogram2d(Lz_30_0, del_30_4, range=[[500,3500], [-500,500]], bins=(80, 60))
    myextent4  =[xedges4[0],xedges4[-1],yedges4[0],yedges4[-1]]
    plt.imshow(H4.T,origin='low',extent=myextent4,interpolation='nearest',aspect='auto', cmap=color)
    
    #messing with the tick labels to make it look pretty
    ax2.set_yticklabels([])
    ax3.set_yticklabels([])
    ax4.set_yticklabels([])
    ax6.set_yticklabels([])
    ax7.set_yticklabels([])
    ax8.set_yticklabels([])
    ax10.set_yticklabels([])
    ax11.set_yticklabels([])
    ax12.set_yticklabels([])
    ax14.set_yticklabels([])
    ax15.set_yticklabels([])
    ax16.set_yticklabels([])
    
    ax1.set_xticklabels([])
    ax5.set_xticklabels([])
    ax9.set_xticklabels([])
    ax4.set_xticklabels([])
    ax8.set_xticklabels([])
    ax12.set_xticklabels([])
    
    ax13.set_xticklabels(['',1000,'',2000,'',3000])
    ax14.set_xticklabels(['',1000,'',2000,'',3000])
    ax15.set_xticklabels(['',1000,'',2000,'',3000])
    ax16.set_xticklabels(['',1000,'',2000,'',3000])
    '''
    ax13.set_xticks([3.0,4.,5.0,6.,7.0,8.,9.0,10.])
    ax14.set_xticks([3.0,4.,5.0,6.,7.0,8.,9.0,10.])
    ax15.set_xticks([3.0,4.,5.0,6.,7.0,8.,9.0,10.])
    ax16.set_xticks([3.0,4.,5.0,6.,7.0,8.,9.0,10.])
    '''
    #adding axis labels for the individual subplots
    ax1.set_ylabel(r'$\theta=15$')
    ax5.set_ylabel(r'$\theta=20$')
    ax9.set_ylabel(r'$\theta=25$')
    ax13.set_ylabel(r'$\theta=30$')
    
    ax13.set_xlabel('t=0.5')
    ax14.set_xlabel('t=1.0')
    ax15.set_xlabel('t=1.5')
    ax16.set_xlabel('t=2.0')
    
    #making vertical lines for CR and lindblad
    ax1.axvline(1760., color='black', alpha=0.5)
    ax1.axvline(2382.25,color='g', ls='dashed', alpha=0.5)
    ax1.axvline(1137.746,color='g', ls='dashed', alpha=0.5)
    ax2.axvline(1760., color='black', alpha=0.5)
    ax2.axvline(2382.25,color='g', ls='dashed', alpha=0.5)
    ax2.axvline(1137.746,color='g', ls='dashed', alpha=0.5)
    ax3.axvline(1760., color='black', alpha=0.5)
    ax3.axvline(2382.25,color='g', ls='dashed', alpha=0.5)
    ax3.axvline(1137.746,color='g', ls='dashed', alpha=0.5)
    ax4.axvline(1760., color='black', alpha=0.5)
    ax4.axvline(2382.25,color='g', ls='dashed', alpha=0.5)
    ax4.axvline(1137.746,color='g', ls='dashed', alpha=0.5)
    ax5.axvline(1760., color='black', alpha=0.5)
    ax5.axvline(2382.25,color='g', ls='dashed', alpha=0.5)
    ax5.axvline(1137.746,color='g', ls='dashed', alpha=0.5)
    ax6.axvline(1760., color='black', alpha=0.5)
    ax6.axvline(2382.25,color='g', ls='dashed', alpha=0.5)
    ax6.axvline(1137.746,color='g', ls='dashed', alpha=0.5)
    ax7.axvline(1760., color='black', alpha=0.5)
    ax7.axvline(2382.25,color='g', ls='dashed', alpha=0.5)
    ax7.axvline(1137.746,color='g', ls='dashed', alpha=0.5)
    ax8.axvline(1760., color='black', alpha=0.5)
    ax8.axvline(2382.25,color='g', ls='dashed', alpha=0.5)
    ax8.axvline(1137.746,color='g', ls='dashed', alpha=0.5)
    ax9.axvline(1760., color='black', alpha=0.5)
    ax9.axvline(2382.25,color='g', ls='dashed', alpha=0.5)
    ax9.axvline(1137.746,color='g', ls='dashed', alpha=0.5)
    ax10.axvline(1760., color='black', alpha=0.5)
    ax10.axvline(2382.25,color='g', ls='dashed', alpha=0.5)
    ax10.axvline(1137.746,color='g', ls='dashed', alpha=0.5)
    ax11.axvline(1760., color='black', alpha=0.5)
    ax11.axvline(2382.25,color='g', ls='dashed', alpha=0.5)
    ax11.axvline(1137.746,color='g', ls='dashed', alpha=0.5)
    ax12.axvline(1760., color='black', alpha=0.5)
    ax12.axvline(2382.25,color='g', ls='dashed', alpha=0.5)
    ax12.axvline(1137.746,color='g', ls='dashed', alpha=0.5)
    ax13.axvline(1760., color='black', alpha=0.5)
    ax13.axvline(2382.25,color='g', ls='dashed', alpha=0.5)
    ax13.axvline(1137.746,color='g', ls='dashed', alpha=0.5)
    ax14.axvline(1760., color='black', alpha=0.5)
    ax14.axvline(2382.25,color='g', ls='dashed', alpha=0.5)
    ax14.axvline(1137.746,color='g', ls='dashed', alpha=0.5)
    ax15.axvline(1760., color='black', alpha=0.5)
    ax15.axvline(2382.25,color='g', ls='dashed', alpha=0.5)
    ax15.axvline(1137.746,color='g', ls='dashed', alpha=0.5)
    ax16.axvline(1760., color='black', alpha=0.5)
    ax16.axvline(2382.25,color='g', ls='dashed', alpha=0.5)
    ax16.axvline(1137.746,color='g', ls='dashed', alpha=0.5)
    
    
    #adding overall axis labels, title, and colorbar label
    fig.text(0.5, 0.02, r'Initial L $(\frac{km^2}{s})$', ha='center', size=18)
    fig.text(0.02, 0.5, r'$\Delta L (\frac{km^2}{s})$', va='center', rotation='vertical', size=18)
    fig.text(0.85, 0.87, 'Number\nof Stars')
    plt.suptitle('Change in Angular Momentum', size=22)
    
    #slight adjustments and put it all together
    fig.subplots_adjust(right=0.8, bottom=0.17, left=0.17)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)
    plt.show()