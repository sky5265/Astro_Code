import numpy as np
import math
from Astro_useful_funcs import *
from Analysis_useful_funcs import *



lat_McDonald = rad_from_degree_min_sec("30:40:18")
print("latitude in rad: "+str(lat_McDonald))

'''
alt = rad_from_degree_min_sec()
print("altitude: "+str(alt))
az = rad_from_degree_min_sec()
print("azimuth: "+str(az))
'''


alts = ["37:45:52", "66:09:40"]
azs = ["50:13:39", "36:07:24"]
time_of_days = ["23:5:56.2", "02:05:56.2"]
days_after_vernal_equinoxs = [88.56, 88.602]


for i in range(len(alts)):
    print(equatorial_from_horizontal(altitude = rad_from_degree_min_sec(alts[i]), azimuth = rad_from_degree_min_sec(azs[i]), latitude = lat_McDonald, time_of_day = time_of_days[i], days_after_vernal_equinox = days_after_vernal_equinoxs[i]))


print("\n\n\n")

print(horizontal_from_equatorial(declination = rad_from_degree_min_sec("48:19:27.7"), RA = radians_from_HHMMSS("19:15:23.63"), latitude = lat_McDonald, time_of_day = "23:05:56.2", days_after_vernal_equinox = 89))
