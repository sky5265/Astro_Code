import sys
sys.path.insert(0, '/Users/sky5265/Documents/GitHub/Astro_Code')
from Astro_useful_funcs import *
from Analysis_useful_funcs import *
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle



def draw_rectangle(x = 0, y = 0, width = 3, height = 2, saveloc = 'rect.png', text = '', rect_color = sKy_colors['blue'], text_color = 'black', fontsize = 30):
    #define Matplotlib figure and axis
    fig, ax = plt.subplots()
    #add rectangle to plot
    ax.add_patch(Rectangle((x, y), width, height, color = rect_color))
    
    x_scale = (width-x)*1.5
    y_scale = x_scale
    
    x_midpoint = x+width/2
    y_midpoint = y+height/2
    
    #plt.xlim([x_midpoint - x_scale/2, x_midpoint + x_scale/2])
    plt.xlim([x, x+width])
    #plt.ylim([y_midpoint - y_scale/2, y_midpoint + y_scale/2])
    plt.ylim([y, y+height])
    
    plt.xticks([])
    plt.yticks([])
    
    ax.spines[['right', 'left', 'top', 'bottom']].set_visible(False)
    plt.text(s = text, x = x_midpoint, y = y_midpoint, horizontalalignment = 'center', verticalalignment = 'center', font = font2, fontsize = fontsize, color = text_color)

    #display plot
    return fig
    
    
    

