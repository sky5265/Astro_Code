import os
import sys
sys.path.insert(0, os.environ['astro_code_dir'])
from Astro_useful_funcs import *
from Analysis_useful_funcs import *

def get_key(dict, value):
    #finds the key for a given value in a dictionary
    for k in dict.keys():
        if dict[k] == value:
            return k

def plot_colors(color_list, saveloc, title = ""):
    min_y = np.inf
    max_y = -np.inf


    plt, _, _ = get_pretty_plot()
    x = np.arange(-10, 10, 0.001)
    for i in range(len(color_list)):
        col = color_list[i]

        b = 5*i + 2
        y = 5*np.sin(0.5*x) + b

        if min(y) < min_y:
            min_y = min(y)
        if max(y) > max_y:
            max_y = max(y)

        plt.plot(x, y, color = col, linewidth = 6)
        plt.text(s = get_key(sKy_colors, col), x = 11, y = y[-1], color = col, fontsize = 15)
    set_pretty_xlabel(plt, "")
    set_pretty_ylabel(plt, "")
    set_pretty_title(plt, title)
    plt.ylim([-5, max_y*1.05])
    plt.xlim([-11, 15])
    plt.xticks([])
    plt.yticks([])
    plt.savefig(saveloc, bbox_inches = 'tight')
    plt.close()



#plotting all colors
plot_colors(color_list = sKy_colors_list, saveloc = 'All Colors.pdf')


#Elsa_frozen
color_list = sKy_color_list(i = 4, scheme="autumn")

plot_colors(color_list = color_list, title = 'BB', saveloc = 'Colors.pdf')









