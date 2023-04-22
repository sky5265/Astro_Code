import sys
sys.path.insert(0, '/Users/sky5265/Documents/GitHub/Astro_Code')
from Astro_useful_funcs import *
from Analysis_useful_funcs import *
import drawing_funcs
from drawing_funcs import *


fig = draw_rectangle(text = 'Engine', rect_color = sKy_colors['mute brown'], text_color = sKy_colors['very dark blue'])
fig.savefig('Engine.pdf', bbox_inches= 'tight')

fig = draw_rectangle(text = 'Shock\nEngine\n\n'+r'$L\sim L_0\left(\frac{t}{t_0}\right)^{-\alpha}$', rect_color = sKy_colors['mute brown'], text_color = sKy_colors['very dark blue'], fontsize = 40)
fig.savefig('Pow_Shock_Engine.pdf', bbox_inches= 'tight')


fig = draw_rectangle(text = 'Radioactive\nEngine\n\n'+r'$m_{Ni, 2}$', rect_color = sKy_colors['mute brown'], text_color = sKy_colors['very dark blue'], width = 3)
fig.savefig('Rad_Decay_Engine_2.pdf', bbox_inches= 'tight')

fig = draw_rectangle(text = 'Radioactive\nEngine\n\n'+r'$m_{Ni, 1}$', rect_color = sKy_colors['mute brown'], text_color = sKy_colors['very dark blue'], width = 3)
fig.savefig('Rad_Decay_Engine_1.pdf', bbox_inches= 'tight')

fig = draw_rectangle(text = 'Radioactive\nEngine\n\n'+r'$m_{Ni}$', rect_color = sKy_colors['mute brown'], text_color = sKy_colors['very dark blue'], width = 3)
fig.savefig('Rad_Decay_Engine.pdf', bbox_inches= 'tight')


fig = draw_rectangle(text = 'Diffusion\n\n'+r'$t_d\sim\sqrt{\frac{m_{ej, 1}}{v}}$', rect_color = sKy_colors['mute brown'], text_color = sKy_colors['very dark blue'])
fig.savefig('Diffusion_1.pdf', bbox_inches= 'tight')

fig = draw_rectangle(text = 'Diffusion\n\n'+r'$t_d\sim\sqrt{\frac{m_{ej, 1}+m_{ej,2}}{v}}$', rect_color = sKy_colors['mute brown'], text_color = sKy_colors['very dark blue'])
fig.savefig('Diffusion_12.pdf', bbox_inches= 'tight')

fig = draw_rectangle(text = 'Diffusion\n\n'+r'$t_d\sim\sqrt{\frac{m_{ej}}{v}}$', rect_color = sKy_colors['mute brown'], text_color = sKy_colors['very dark blue'])
fig.savefig('Diffusion.pdf', bbox_inches= 'tight')


fig = draw_rectangle(text = 'Blackbody\nSED+\nPhotometry', rect_color = sKy_colors['mute brown'], text_color = sKy_colors['very dark blue'])
fig.savefig('SED.pdf', bbox_inches= 'tight')


fig = draw_rectangle(text = 'Blackbody', rect_color = sKy_colors['mute brown'], text_color = sKy_colors['very dark blue'], fontsize = 40)
fig.savefig('Blackbody.pdf', bbox_inches= 'tight')

fig = draw_rectangle(text = 'Photometry', rect_color = sKy_colors['mute brown'], text_color = sKy_colors['very dark blue'], fontsize = 40)
fig.savefig('Photometry.pdf', bbox_inches= 'tight')
