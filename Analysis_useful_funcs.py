import numpy as np
import math
import matplotlib.pyplot as plt

def find_nearest(wv,value):
    '''inputs: wv-list of data to look in, value-the value we are looking for in wv
       outputs: idx-index in wv that is closest to value
       this finds the index of the closest element to value in list wv'''
    idx= (np.abs(np.asarray(wv) - value)).argmin()
    return idx

def pretty_plot(x, y, xlabel = r'', ylabel = '', title = '', label = '', xlim = [], ylim = [], save_loc='', display_or_nah = False):
    '''inputs: x-indpendent variable of plot to create, y-dependent variable of plot to create, xlabel-string of label on x-axis, ylabel-string of label on y-axis, title-string of title for figure, xlim-list of limits of plot for x-axis, ylim-list of limits of plot for y-axis, labels-string of label of plot, save_loc-string location in file where to store the plot (default to ''; if '' is set, the figure won't be saved), display_or_nah-boolean variable to control whether the figure should be displayed (default to False)
        outputs: none
        This function plots x and y using plt.plot--I only have this function because I have a specific way I like plotting, and hate remembering all the details, so this function plots things for me'''
    W = 15
    H = 10
    plt.figure(figsize=(W, H))
    
    if len(label) > 0:
        plt.plot(x, y, linewidth=H/3, label=label)
    else:
        plt.plot(x, y, linewidth=H/3)
    if len(xlabel) > 0:
        plt.xlabel(xlabel, fontsize=2*W)
        
    plt.xticks(fontsize=1.5*W)
    
    if len(xlim) > 0:
        plt.xlim(xlim)
        
    if len(ylim) > 0:
        plt.ylim(ylim)
    
    if len(ylabel) > 0:
        plt.ylabel(ylim, fontsize=2.5*H)
        
    plt.yticks(fontsize=1.5*W)
    
    if len(title) > 0:
        plt.title(title, fontsize=2*W)
        
    if len(label):
        plt.legend(fontsize = 1.5*H)
        
    if len(save_loc) > 0:
        plt.savefig(save_loc)
        
    if display_or_nah:
        plt.show()
