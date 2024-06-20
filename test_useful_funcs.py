import os
import sys
sys.path.insert(0, os.environ['astro_code_dir'])
from Astro_useful_funcs import *
from Analysis_useful_funcs import *

def plot_colors(color_list, saveloc, keys = [], title = ""):
    min_y = np.inf
    max_y = -np.inf


    plt, _, _ = get_pretty_plot()
    x = np.arange(-10, 10, 0.01)

    offsets = np.random.random(len(x))*3

    for i in range(len(color_list)):
        col = color_list[i]

        b = 2*i + 1
        y = (5*np.sin(0.5*x) + b)

        if min(y) < min_y:
            min_y = min(y)
        if max(y) > max_y:
            max_y = max(y)

        #plt.scatter(x, y+offsets, color = col, s = 100, edgecolor = sKy_colors['grey'])
        plt.scatter(x, y+offsets, color = col, s = 100)
        #plt.plot(x, y, color = col, linewidth = 6, mec = sKy_colors['grey'], mew = 3)
        if len(keys) > i:
            plt.text(s = str(keys[i]), x = 11, y = y[-1], color = col, fontsize = 15)
    set_pretty_xlabel(plt, "")
    set_pretty_ylabel(plt, "")
    set_pretty_title(plt, title)
    plt.ylim([-5, max_y*1.05])
    plt.xlim([-11, 15])
    #plt.plot([-99999, 9999], [0, 0], linewidth = 100000, color = 'black', zorder = -10)
    plt.xticks([])
    plt.yticks([])
    plt.savefig(saveloc, bbox_inches = 'tight')
    plt.close()



#plotting all colors

plot_colors(color_list = [(1,0,0), (0,1,0), (0,0,1)], title = '', 
    keys = [str((1,0,0)), str((0,1,0)), str((0,0,1))], saveloc = '3 Colors.pdf')




mkdir("Schemes/")
scheme_name = []
color_list = []
for scheme in color_schemes:
    

    plot_colors(color_list = color_schemes[scheme], title = scheme, saveloc = "Schemes/"+scheme+".pdf")

mkdir("New Schemes/")

sonmi = [(0.9, 0.1, 0.7), (0.6, 0.1, 0.7), (0.3, 0.1, 0.7)]
plot_colors(color_list = sonmi, title = '', 
    keys = sonmi, saveloc = 'New Schemes/Sonmi.pdf')

def exaggerate_color(color, factor = 2.0):
    return (min(color[0]**(factor), 1.0), min(color[1]**factor, 1.0), min(color[2]**factor, 1))

def exaggerate_green(color, factor = 2.0):
    return (color[0], min(color[1]**factor, 1.0), color[2])

def transformation(colors, function, factor = 2.0):
    tor = []
    for col in colors:
        tor.append(function(col, factor = factor))
    return tor

plot_colors(color_list = transformation(colors = color_schemes['rainforest_flower'], 
    function = exaggerate_color, factor = 2.0), 
    title = 'Exaggerated rainforest_flower', saveloc = 'Exaggerated Color.pdf')

plot_colors(color_list = transformation(colors = color_schemes['rainforest_flower'], 
    function = exaggerate_color, factor = 0.5), 
    title = 'Muted rainforest_flower', saveloc = 'Muted Color.pdf')


color_list = transformation(colors = color_schemes['rainforest_flower'], 
    function = exaggerate_green, factor = 2.0)
plot_colors(color_list = color_list, 
    title = 'Muted Green rainforest_flower', saveloc = 'Muted Green Color.pdf')


manitoba = [
(0.8, 0.9, 0.1), (0.8, 0.6, 0.1), 
(0.2, 0.8, 0.1), sKy_colors['grey'], sKy_colors['orange']
]

#how to add a black border around the dots in scatter?





'''

def get_rand(min, max):
    num = random.random()
    num *= (max-min) + min
    return num

def color_generator(red_window, green_window, blue_window):
    red = get_rand(red_window[0], red_window[1])
    green = get_rand(green_window[0], green_window[1])
    blue = get_rand(blue_window[0], blue_window[1])

    return (red, green, blue)


def scheme_generator(num_colors, red_window, green_window, blue_window):

    colors = []
    keys = []
    for i in range(num_colors):
        color = color_generator(red_window, green_window, blue_window)
        color = (round(color[0], 2), round(color[1], 2), round(color[2], 2))
        colors.append(color)
        keys.append(str(color))

    return colors, keys



red_window = [0.8, 1.0]
green_window = [0.4, 0.6]
blue_window = [0, 0.1]
colors, keys = scheme_generator(num_colors = 10, red_window = red_window, 
    green_window = green_window, blue_window = blue_window)

plot_colors(color_list = colors, keys = keys, title = 'New Scheme', saveloc = 'New Schemes/1.pdf')
'''




'''


mkdir("colors/")

reds = np.arange(0.2, 1.1, 0.2)
greens = np.arange(0.2, 1.1, 0.2)
blues = np.arange(0.2, 1.1, 0.2)


for red in reds:
    color_list = []
    keys = []
    for green in greens:
        color = (round(red, 5), round(green, 5), 0.1)
        color_list.append(color)
        keys.append(str(color))

    for green in greens:
        color = (round(red, 5), round(green, 5), 0.5)
        color_list.append(color)
        keys.append(str(color))

    for green in greens:
        color = (round(red, 5), round(green, 5), 0.8)
        color_list.append(color)
        keys.append(str(color))


    plot_colors(color_list = color_list, title = '', keys = keys, saveloc = 'colors/Greens (R = '+str(round(red, 2))+').pdf')

for green in greens:
    color_list = []
    keys = []
    for red in reds:
        color = (round(red, 5), round(green, 5), 0.1)
        color_list.append(color)
        keys.append(str(color))

    for red in reds:
        color = (round(red, 5), round(green, 5), 0.5)
        color_list.append(color)
        keys.append(str(color))

    for red in reds:
        color = (round(red, 5), round(green, 5), 0.8)
        color_list.append(color)
        keys.append(str(color))


    plot_colors(color_list = color_list, title = '', keys = keys, saveloc = 'colors/Reds (G = '+str(round(green, 2))+').pdf')

for blue in blues:
    color_list = []
    keys = []
    for red in reds:
        color = (round(red, 5), 0.1, round(blue, 5))
        color_list.append(color)
        keys.append(str(color))

    for red in reds:

        color = (round(red, 5), 0.5, round(blue, 5))
        color_list.append(color)
        keys.append(str(color))

    for red in reds:

        color = (round(red, 5), 0.8, round(blue, 5))
        color_list.append(color)
        keys.append(str(color))



    plot_colors(color_list = color_list, title = '', keys = keys, saveloc = 'colors/Reds (B = '+str(round(blue, 2))+').pdf')


'''
