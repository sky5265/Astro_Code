import numpy as np
import math
import matplotlib.pyplot as plt
import os


def isint(element):
    return is_int(element)
    
def is_int(element):
    '''checks if an element can be turned into an int'''
    try:
        int(element)
        return True
    except ValueError:
        return False

def isfloat(element):
    return is_float(element)


def write_to_file(file_loc, data_out, separator = ' ', headers=None, append=False):
    '''inputs: file_loc-location to which to write to, data_to_write-this is a 2d array that will be written to a file
       outputs: none
       This function writes a 2d array to a file'''
       
    if not append:
        out = open(file_loc ,'w')
    else:
        out = open(file_loc ,'a+')
    
    if type(data_out) is str:
        out.write(data_out+"\n")
        
    else:
        for i  in range(len(data_out)):
            if isinstance(data_out[i], float) or type(data_out[i]) is str:
                out.write(str(data_out[i])+separator)
            else:
                for j in range(len(data_out[i])):
                    out.write(str(data_out[i][j])+separator)
                
            out.write("\n")
    out.close()
            
            

def load_file(file_loc, header_marker=None):
    '''same as read_file'''
    
    return read_file(file_loc, header_marker = header_marker)

def read_file(file_loc, comments ='#', separators = [], skip_lines = [-1], header_row = None, header_marker=None):
    '''inputs: file_loc-string of file location to read in
       outputs: data-2d array of data read in and passed out
       This function reads in a dataset in a file located at file_loc and returns a 2d array of this'''
    headers = ''
    row_num = 0
    it = 0
    skipped_lines = []
    with open(file_loc) as f:
        lines = []
        for line in f:
            if row_num not in skip_lines:
                if len(separators) > 0:
                    for separator in separators:
                        line = line.replace(separator, ' ')
                if len(line.split()) > 0 and line.split()[0] != comments:
                    if header_row is not None and it == 0:
                        headers = line.rstrip()
                        
                    else:
                        lines.append(line.rstrip().split())
                    it += 1
            elif row_num in skip_lines:
                skipped_lines.append(line)
            row_num += 1
        
    if len(separators) > 0:
        for separator in separators:
            headers = headers.replace(separator, ' ')
    
    if len(headers) > 0:
        headers = headers.split()
    
    dict = {}
    if len(headers) > 0:
        for header in headers:
            dict[header] = []
    
    for i in range(len(lines)):
        for j in range(len(lines[i])):
            if len(headers) < j+1:
                if j not in dict.keys():
                    dict[j] = []
                dict[j].append(lines[i][j])
            else:
                if headers[j] not in dict.keys():
                    dict[headers[j]] == []
                dict[headers[j]].append(lines[i][j])
        
    return {"lines": lines, "dict": dict, "skipped_lines":skipped_lines}
        
    

def is_float(element):
    '''checks if an element can be turned into a float'''
    try:
        float(element)
        return True
    except ValueError:
        return False


def find_nearest(wv,value):
    '''inputs: wv-list of data to look in, value-the value we are looking for in wv
       outputs: idx-index in wv that is closest to value
       this finds the index of the closest element to value in list wv'''
    idx= (np.abs(np.asarray(wv) - value)).argmin()
    return idx

def mkdir(dir_loc):
    '''inputs: dir_loc-location of directory to create
       outputs: none
       This function creates the directory at dir_loc if it doesn't exist. If it does exist, nothing happens'''
    if not os.path.isdir(dir_loc):
        os.mkdir(dir_loc)
    
    

def pretty_plot(x, y, xlabel = r'', ylabel = '', title = '', label = '', xlim = [], ylim = [], xscale='', yscale='', save_loc='', display_or_nah = False):
    '''inputs: x-indpendent variable of plot to create, y-dependent variable of plot to create, xlabel-string of label on x-axis, ylabel-string of label on y-axis, title-string of title for figure, xlim-list of limits of plot for x-axis, ylim-list of limits of plot for y-axis, labels-string of label of plot, save_loc-string location in file where to store the plot (default to ''; if '' is set, the figure won't be saved), display_or_nah-boolean variable to control whether the figure should be displayed (default to False)
        outputs: none
        This function plots x and y using plt.plot--I only have this function because I have a specific way I like plotting, and hate remembering all the details, so this function plots things for me'''
    W = 15
    H = 10

    
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
        plt.ylabel(ylabel, fontsize=2.5*H)
        
    plt.yticks(fontsize=1.5*W)
    
    if len(title) > 0:
        plt.title(title, fontsize=2*W)
        
    if len(label):
        plt.legend(fontsize = 1.5*H)
    
    if xscale == 'log':
        plt.xscale('log')
        
    if yscale == 'log':
        plt.yscale('log')
        
    if len(save_loc) > 0:
        plt.savefig(save_loc)
        
    if display_or_nah:
        plt.show()
        
    return plt
