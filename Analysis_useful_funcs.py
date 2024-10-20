import numpy as np
import math
import matplotlib.pyplot as plt
import os
import matplotlib
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import shutil
import pandas as pd
import sys




sKy_colors = {'light blue':'#63B8FF', 'blue':'#4876FF', 'very dark blue':'#27408B', 
'blue grey':'#C6E2FF', 'dim cyan':'#98F5FF', 'cyan':'#00FFFF','red':'#FF4040', 
'mute red':'#EE6363', 'dark mute red':'#CD5555', 'dark red':'#CD2626', 'green':'#00FF7F', 
'honest green':'#008B45', 'dark green':'#008B45', 'grey':'#8B8989', 'dark grey':'#666666', 
'orange':'#FF9912', 'purple':'#8E388E', 'magenta':'#FF00FF', 'purple pink':'#FF83FA', 
'dark purple pink':'#BF3EFF', 'bright brown':'#8B5A00', 'dull brown':'#8B4726', 'mute brown':'#BC8F8F',
'light grey':(0.8, 0.8, 0.8)}

color_schemes = {
"elsa":[sKy_colors["blue"], sKy_colors["blue grey"], 
sKy_colors["dark purple pink"], sKy_colors["cyan"], sKy_colors["light blue"], 
sKy_colors["magenta"], sKy_colors["grey"]],

"simplex":[sKy_colors['orange'], sKy_colors['grey'], sKy_colors['mute brown'], sKy_colors['blue grey']],

"autumn": [sKy_colors['grey'], sKy_colors['mute red'], sKy_colors['orange'], sKy_colors['dark red']],

"2145": [sKy_colors['dark purple pink'], 
sKy_colors['light blue'], 
sKy_colors['purple pink'],  sKy_colors['blue'],
 sKy_colors['magenta'],  
 sKy_colors['very dark blue']],

"ocean": [(0.8, 0.5, 1.0), (0.8, 0.8, 1.0), 
(0.2, 0.5, 1.0), (0.3, 0.0, 1.0),
(0.2, 0.8, 1.0),  (0.2, 1.0, 0.8),
(0.6, 0.8, 0.6)],

"red_blue": [sKy_colors['dark red'], 
sKy_colors['mute red'], (1.0, 0.6, 0.8),
sKy_colors['orange'], 
sKy_colors['very dark blue'], 
sKy_colors['blue'], sKy_colors['light blue'],
(0.2, 0.1, 0.6)],

"rainforest_flower":[(0.2, 0.6, 0.5),
(0.8, 0.8, 0.8), (1.0, 0.4, 0.3),( 1.0, 0.5, 0.8),
(1.0, 0.8, 0.1), (1.0, 0.2, 0.1) ],

"childhood":[sKy_colors['dark purple pink'], 
sKy_colors['orange'], sKy_colors['magenta'], 
 sKy_colors['green'], sKy_colors['blue'], 
 sKy_colors['red'], sKy_colors['cyan'],
 (1.0, 1.0, 0.5)],

"chill": [sKy_colors['light blue'], 
sKy_colors['blue grey'], sKy_colors['mute red'], 
sKy_colors['grey'], sKy_colors['purple pink'], 
sKy_colors['mute brown'], (0.8, 1.0, 0.8),
(1.0, 0.6, 0.5)],

"icecream": ['black', 'pink', sKy_colors['grey'],  sKy_colors['purple pink']],

"0605":[(0.1, 0.6, 0.5), 
(0.7, 0.6, 0.5), (0.9, 0.6, 0.5), 
sKy_colors['grey'], (1.0, 0.8, 0.1), 
(1.0, 0.6, 0.8)
],

"space_person":[(0.26099, 0.55795, 0.98008), 
(0.44243, 0.91793, 0.85), 
(0.49594, 0.15375, 0.5658), 
(0.3246, 0.20608, 0.53124), 
(0.14722, 0.43899, 0.49572), 
(0.21402, 0.19033, 0.45783)]

}


sKy_colors_list = [sKy_colors[i] for i in sKy_colors.keys()]

font1 = 'Shree Devanagari 714' #use this for unbolded text
font2 = 'hiragino sans' #use this for bolded text
W = 15
H = 10
fig, ax = plt.subplots(figsize=(W, H))

#ax.xaxis.set_major_locator(MultipleLocator(20))
#ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))

# For the minor ticks, use no labels; default NullFormatter.
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
plt.tick_params(
    which='major',
    bottom='on',
    top='on',
    left='on',
    right='on',
    direction='in',
    length=25)
plt.tick_params(
    which='minor',
    bottom='on',
    top='on',
    left='on',
    right='on',
    direction='in',
    length=10)
matplotlib.rcParams["font.family"] = 'serif'
#matplotlib.rc('axes', unicode_minus=False)
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'serif'
matplotlib.rcParams['mathtext.it'] = 'serif'
matplotlib.rcParams['mathtext.bf'] = 'serif'
plt.xticks(fontsize=1.5*W)
plt.yticks(fontsize=1.5*W)

def get_element_symbol(atomic_number):
    df = pd.read_csv(os.environ['astro_code_dir']+'/Periodic_Table.csv')
    element = df.loc[df['AtomicNumber'] == int(atomic_number)]
    l = list(element['Symbol'])
    if len(l) == 1:
        return l[0]
    else:
        raise ValueError("Atomic Number: "+str(atomic_number)+" is invalid!")

def rm(file):
    '''inputs: file-string file name of file to be removed
       outputs: none
       This function deletes a file or folder given the string of the location'''
    if os.path.isdir(file):
        os.system("rm -r "+file)
    elif os.path.isfile(file):
        os.remove(file)


def extract_digits(filename):
    # Use regex to find all digits in the filename
    digits = re.findall(r'\d+', filename)
    # Join the digits and convert to an integer
    return int(''.join(digits)) if digits else None


def get_key(dict, value):
    #finds the key for a given value in a dictionary
    for k in dict.keys():
        if dict[k] == value:
            return k

def get_colors(i, scheme = None):
    return sKy_color_list(i, scheme)

def sKy_color_list(i, scheme = None):
    '''inputs: i-integer number of colors to generate, scheme-name of color scheme to choose colors from
       outputs: color_list-list of colors generated from sKy_colors dictionary with length i
       This function generates a list of i colors that can be used in a color map. It doesn't 
       do anything super intelligent. It just picks i random colors from the dictionary. Scheme allows you to 
       specify if you want the colors to come from a color scheme that exists. It will pick the first i colors 
       from the specified color scheme, if the scheme exists. If i is larger than the number of colors in
       the scheme, more colors will be chosen from the rest of the set of colors
    '''

    import random
    tor = None
    if scheme is None or scheme not in color_schemes.keys():
        idx_list = np.asarray(random.sample(range(len(sKy_colors.keys())), i))
        tor = [sKy_colors_list[i] for i in idx_list]
    else:

        tor = [color_schemes[scheme][j] for j in range(min([i, len(color_schemes[scheme])]))]
        if len(tor) < i:
            idx_list = np.asarray(random.sample(range(len(sKy_colors.keys())), i-len(tor)))
            for k in idx_list:
                tor.append(sKy_colors_list[k])

    return tor


def plot_available_fonts(save_loc = 'fonts.png', bold = False):
    import matplotlib.font_manager
    fpaths = matplotlib.font_manager.findSystemFonts()

    all_fonts = []

    for i in fpaths:
        
        try:
            f = matplotlib.font_manager.get_font(i)
            print(f.family_name)
            all_fonts.append(f.family_name)
        except:
            pass
    W = 300
    H = 20
    xs_added = []
    values_to_skip = []
    #values_to_skip = [1, 2, 4, 5, 7, 8, 9, 11, 12, 13, 14, 15, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 29, 30, 31, 32, 37, 41, 43, 44, 48, 49, 50, 52, 55, 56, 59, 62, 65, 67, ]
    plt.figure(figsize=(W, H))
    for i in range(len(all_fonts)):
        if i in values_to_skip:
            continue
        try:
            xs_added.append((math.floor(i*0.1))*160.0)
            if bold:
                plt.text(x = (math.floor(i*0.1))*160.0, y = (i*0.1)%1, s = all_fonts[i]+": Absolute Magnitude", font_properties = all_fonts[i], fontsize = 20, weight = 'bold')
            else:
                plt.text(x = (math.floor(i*0.1))*160.0, y = (i*0.1)%1, s = all_fonts[i]+": Absolute Magnitude", font_properties = all_fonts[i], fontsize = 20)
            plt.text(x = (math.floor(i*0.1))*160.0, y = (i*0.1)%1+0.03, s = str(i), fontsize = 15)
        except:
            pass

    plt.xlim([0, 10000])
    plt.savefig(save_loc)

def isint(element):
    return is_int(element)
    
def file_exists(file):
    return os.path.exists(file)
    
def is_int(element):
    '''checks if an element can be turned into an int'''
    try:
        int(element)
        return True
    except ValueError:
        return False

def isfloat(element):
    return is_float(element)
    
    
def conv_to_deluxe_table(file_loc, deluxe_table_loc = None, table_caption = None, table_label = None, separators = [' ', ','], insert_columns = None, switch_columns = None, insert_rows = None):
    '''inputs: file_loc-location of file with data to turn into latex code for deluxe table, deluce_table_loc-location of where to write code for deluxe table to, insert_columns-l-the number of columns to be inserted-this will create this many columns near the end of the table, with randomized column names, switch_columns-will move columns around-this should be a list like["1,2"]-which will tell it to switch columns 1 and 2-indexing from 0, insert_rows-number of rows to be inserted-this will be done after all the column switching and adding is done
       outputs: none
       this will take a location of a file that holds data and convert it into latex code for a deluxe table. If deluxe_table_loc is passed, the code will be written to that location
    '''
    
    random_col_names = ['Cheese', 'Lisa', 'Jennie', 'Pizza', 'Bibimbap', 'PaniPuri', 'PavBhaji', 'Fries', 'Dip', 'Hummus', 'Pita']
    data = read_file(file_loc, separators = separators)['dict']
    
    keys = list(data.keys())
    nrow = len(data[keys[0]])
    if insert_columns is not None:
        for i in range(insert_columns):
            if i+len(keys) < len(random_col_names)-1:
                colname = random_col_names[i+len(keys)]
            else:
                colname = "Column "+str(i)
            data[colname] = [str(i) for i in np.zeros(nrow)]
            
    keys = list(data.keys())
    
    if switch_columns is None:
        switch_columns = []
    for switch in switch_columns:
        a = keys[int(switch[:switch.index(',')])]
        b = keys[int(switch[switch.index(',')+1:])]
        
        temp = [i for i in data[a]]
        data[a] = [i for i in data[b]]
        data[b] = [i for i in temp]
    
    keys = list(data.keys())
    ncol = len(keys)
    
    latex_code = "\\begin{deluxetable}{"+"c"*ncol+"}\n"
    latex_code += "\\tabletypesize{\\footnotesize}\n"
    latex_code += "\\tablecolumns{"+str(ncol)+"}\n"
    
    if table_caption is None:
        latex_code += "\\tablecaption{I like Cheese. "
    else:
        latex_code += "\\tablecaption{"+str(table_caption)+". "
    if table_label is None:
        latex_code += "\\label{table:fun_table}.}\n"
    else:
        latex_code += "\\label{table:"+str(table_label)+"}.}\n"
        
    col_heads = ""
    for i in range(ncol-1):
        if i > len(random_col_names)-1:
            col_name = keys[i]
        else:
            col_name = random_col_names[i]
        col_heads += " \\colhead{"+str(col_name)+"} &"
        
    if i > len(random_col_names)-1:
        col_name = keys[i]
    else:
        col_name = random_col_names[i]
    
    col_heads += " \\colhead{"+str(col_name)+"} "
    
    latex_code += "\\tablehead{"+str(col_heads)+"}\n"
    latex_code += "\startdata\n"
    
    for i in range(nrow):
        #read the value at each row from the keys
        line = ''
        for j in range(len(keys)-1):
            key = keys[j]
            line += data[key][i] + "\t&\t"
        #print("keys: "+str(keys))
        #print("keys[len(keys)-1]: "+str(keys[len(keys)-1]))
        line += data[keys[-1]][i]+" \\\\\n"
        latex_code += line
        
    latex_code += '\enddata\n'
    latex_code += '\end{deluxetable}'
    
    if deluxe_table_loc is not None:
        write_to_file(deluxe_table_loc, latex_code)
    else:
        print(latex_code)
    
    
    
    
    
def ls(dir = None, extension = None):
    '''inputs: dir-directory to read the list of files from, extension-which extension should the files have to be listed (default to none-will list all file types)
       outputs: list of files WITHOUT directory name
       This function gets a list of all files within a directory name. Optionally, can input a specific file extension to only list out files with that extension. This will only list out the short name of the files, like ls would on terminal
    
    '''
    
    files_list = []
    if dir is not None:
    
    
        for filename in os.listdir(dir):
            ext = filename[-1*len(extension):]
            short_name = filename
            if extension is not None and ext == extension:
                files_list.append(short_name)
    else:
    
        for filename in os.listdir():
            ext = filename[-1*len(extension):]
            short_name = filename
            if extension is not None and ext == extension:
                files_list.append(short_name)
        
    return files_list
            


def write_to_file(file_loc, data_out, separator = ' ', headers=None, append=False):
    '''inputs: file_loc-location to which to write to, data_to_write-this is a 2d array that will be written to a file
       outputs: none
       This function writes a 2d array to a file'''
       
    if not append:
        out = open(file_loc ,'w')
    else:
        out = open(file_loc ,'a+')
    
    if type(data_out) is str:
        out.write(data_out)
        
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

def read_file(file_loc, comments ='#', separators = [], skip_lines = [-1], header_list = None, header_row = None, header_line = None, header_marker=None, find_floats = False):
    '''inputs: file_loc-string of file location to read in, comments-which character makes the line a comment, separators-list of characters that are used to separate line entry from line entry, skip_lines-list of lines to skip when reading in, header_row-which row in the datafile is the set of headers (defaults to none, if none, the dictionary will have keys of just row numbers), header_line-aame as header_row, header_marker-the character that marks a line as the header, find_floats-will only return a data entry in the final dictionary if it is a float
       outputs: data-2d array of data read in and passed out
       This function reads in a dataset in a file located at file_loc and returns a dictionary of this'''
       
    if header_line != None and header_row == None:
        header_row = header_line
        
        
    
        
    headers = ''
    if header_list is not None:
        headers = header_list.rstrip()
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
                    if len(headers) == 0 and header_row is not None and it == 0:
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
    except:
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
def cp(src, dst):
    '''inputs: src-name of source file to copy from
       outputs: dst-name of destination file to copy to
       This function copies a file from src to dst'''
    shutil.copy(src, dst)


def mv(src, dst):
    '''inputs: src-name of source file to copy from
       outputs: dst-name of destination file to copy to
       This function copies a file from src to dst'''
    if os.path.isfile(src) and (os.path.isfile(dst) or os.path.isdir(dst)):
        shutil.move(src, dst)
        
def get_order_of_mag(num):
    '''returns order of magnitude for a number passed in. This returns just the exponent, so if you pass in 8.9E-45, it will return -45'''
    
    if num == 0:
        return 0

    absnum = abs(num)
    order = math.log10(absnum)
    res = math.floor(order)

    return res


def set_pretty_title(plt, title, fontsize = 2*W, weight = 'normal'):
    plt.title(title, fontsize=fontsize, weight = weight)
    return plt

def set_pretty_xlabel(plt, xlabel, fontsize = 2*W, weight = 'normal'):
    plt.xlabel(xlabel, fontsize=fontsize, weight = weight)
    return plt

def set_pretty_ylabel(plt, ylabel, fontsize = 2*W, weight = 'normal'):
    plt.ylabel(ylabel, fontsize=fontsize, weight = weight)
    return plt

def set_pretty_legend(plt, fontsize = 1.5*H, ncol = np.nan):
    if ncol == np.nan:
        plt.legend(fontsize = fontsize)   
    else:
        plt.legend(fontsize = fontsize, ncol = ncol)
    return plt


def get_pretty_plot(W = 15, H = 10):
    fig, ax = plt.subplots(figsize=(W, H))

    #ax.xaxis.set_major_locator(MultipleLocator(20))
    #ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))

    # For the minor ticks, use no labels; default NullFormatter.
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.tick_params(
        which='major',
        bottom='on',
        top='on',
        left='on',
        right='on',
        direction='in',
        length=25)
    plt.tick_params(
        which='minor',
        bottom='on',
        top='on',
        left='on',
        right='on',
        direction='in',
        length=10)
    matplotlib.rcParams["font.family"] = 'serif'
    #matplotlib.rc('axes', unicode_minus=False)
    matplotlib.rcParams['mathtext.fontset'] = 'custom'
    matplotlib.rcParams['mathtext.rm'] = 'serif'
    matplotlib.rcParams['mathtext.it'] = 'serif'
    matplotlib.rcParams['mathtext.bf'] = 'serif'
    plt.xticks(fontsize=1.5*W)
    plt.yticks(fontsize=1.5*W)
    weight = 'normal'
    
    '''
    if dark_theme:
        plt.style.use('dark_background')
    elif minimalist:
        plt.style.use('fivethirtyeight')
    else:
        plt.style.use('seaborn')'''
    return plt, W, H

    

def pretty_plot(x, y, xlabel = r'', ylabel = '', title = '', label = '', xlim = [], ylim = [], xscale='', yscale='', vlines = [], hlines = [], save_loc='', display_or_nah = False, scatter = False, s = None, dark_theme = False, minimalist = False, color = None, bold = False):
    '''inputs: x-indpendent variable of plot to create, y-dependent variable of plot to create, xlabel-string of label on x-axis, ylabel-string of label on y-axis, title-string of title for figure, xlim-list of limits of plot for x-axis, ylim-list of limits of plot for y-axis, labels-string of label of plot, save_loc-string location in file where to store the plot (default to ''; if '' is set, the figure won't be saved), display_or_nah-boolean variable to control whether the figure should be displayed (default to False)
        outputs: none
        This function plots x and y using plt.plot--I only have this function because I have a specific way I like plotting, and hate remembering all the details, so this function plots things for me'''
    
        
    W = 15
    H = 10
    fig, ax = plt.subplots(figsize=(W, H))

    #ax.xaxis.set_major_locator(MultipleLocator(20))
    #ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))

    # For the minor ticks, use no labels; default NullFormatter.
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.tick_params(
        which='major',
        bottom='on',
        top='on',
        left='on',
        right='on',
        direction='in',
        length=25)
    plt.tick_params(
        which='minor',
        bottom='on',
        top='on',
        left='on',
        right='on',
        direction='in',
        length=10)
    matplotlib.rcParams["font.family"] = 'serif'
    #matplotlib.rc('axes', unicode_minus=False)
    matplotlib.rcParams['mathtext.fontset'] = 'custom'
    matplotlib.rcParams['mathtext.rm'] = 'serif'
    matplotlib.rcParams['mathtext.it'] = 'serif'
    matplotlib.rcParams['mathtext.bf'] = 'serif'

    weight = 'normal'
    
    '''
    if dark_theme:
        plt.style.use('dark_background')
    elif minimalist:
        plt.style.use('fivethirtyeight')
    else:
        plt.style.use('seaborn')'''
        
    
    if color is None:
        color = sKy_colors[list(sKy_colors.keys())[0]]
    
    if not scatter:
    
        if len(label) > 0:
            plt.plot(x, y, linewidth=H/3, label=label, color = color)
        else:
            plt.plot(x, y, linewidth=H/3, color = color)
    else:
        if s == None:
            s = 5
        if len(label) > 0:
            plt.scatter(x, y, s = H*s, label=label, color = color)
        else:
            plt.scatter(x, y, s = H*s, color = color)



    if len(xlabel) > 0:
        plt.xlabel(xlabel, fontsize=2*W, weight = weight)
        
    plt.xticks(fontsize=1.5*W)
    
    if len(xlim) > 0:
        plt.xlim(xlim)
        
    if len(ylim) > 0:
        plt.ylim(ylim)
    
    if len(ylabel) > 0:
        plt.ylabel(ylabel, fontsize=2*W, weight = weight)
        
    plt.yticks(fontsize=1.5*W)

    if len(vlines) > 0:
        ylim = plt.gca().get_ylim()
        plt.vlines(vlines,ylim[0], ylim[1])

    if len(hlines) > 0:
        xlim = plt.gca().get_xlim()
        plt.hlines(hlines,xlim[0], xlim[1])

    if len(title) > 0:
        set_pretty_title(plt = plt, title = title, weight = weight, fontsize = 2*W)
        
    if len(label):
        plt.legend(fontsize = 1.5*H, weight = weight)
    
    if xscale == 'log':
        plt.xscale('log')
        
    if yscale == 'log':
        plt.yscale('log')
        
    if len(save_loc) > 0:
        plt.savefig(save_loc, bbox_inches=  'tight')
        
    if display_or_nah:
        plt.show()
        
    plt.close()
        
    return plt
