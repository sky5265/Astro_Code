import numpy as np
import math
from Astro_useful_funcs import *
from Analysis_useful_funcs import *



galactic_l_deg = 144.5
galactic_b_deg = 6.5

old_eq_coord = equatorial_from_galactic(b_deg = galactic_b_deg, l_deg = galactic_l_deg)




curr_coord = precess_ecliptic_coord_fwd(num_years_to_precess = 72, lam_old_deg = galactic_l_deg, beta_old_deg = galactic_b_deg)
curr_lam_deg = curr_coord['lam_new_deg']
curr_beta_deg = curr_coord['beta_new_deg']

eq_coord = equatorial_from_galactic(b_deg = curr_beta_deg, l_deg = curr_lam_deg)

print('eq coord: '+str(eq_coord))






lat_McDonald_rad = rad_from_degree_min_sec("30:40:18")
lat_McDonald_deg = lat_McDonald_rad*180./(np.pi)

long_McDonald_rad = radians_from_HHMMSS("-6:56:5.2")
long_McDonald_deg = long_McDonald_rad*180./(np.pi)

lat_AA_rad = rad_from_degree_min_sec("-31:16:37")
lat_AA_deg = lat_AA_rad*180./np.pi

long_AA_rad = radians_from_HHMMSS("9:56:15.9")
long_AA_deg = long_AA_rad*180./np.pi

print("latitude in deg: "+str(lat_McDonald_deg))
print("longitude in deg: "+str(long_McDonald_deg))


'''
alt = rad_from_degree_min_sec()
print("altitude: "+str(alt))
az = rad_from_degree_min_sec()
print("azimuth: "+str(az))
'''


alts = ["-16:34:30", "-8:32:10"]
azs = ["138:57:20", "144:23:34"]
time_of_days = ["12:57:25.5", "13:57:25.5"]
days_after_vernal_equinoxs = [89, 89]

'''
alts = ["8:18:01", "4:17:40"]
azs = ["346:29:54", "337:18:12"]
time_of_days = ["2:57:25.5", "3:57:25.5"]
days_after_vernal_equinoxs = [90, 90]
'''

for i in range(len(alts)):
    print(equatorial_from_horizontal(altitude = rad_from_degree_min_sec(alts[i]), azimuth = rad_from_degree_min_sec(azs[i]), latitude_deg = lat_McDonald_deg, time_of_day = time_of_days[i], days_after_vernal_equinox = days_after_vernal_equinoxs[i], longitude_deg = long_McDonald_deg, UTC_offset = -5))
    
    
    #print(equatorial_from_horizontal(altitude = rad_from_degree_min_sec(alts[i]), azimuth = rad_from_degree_min_sec(azs[i]), latitude_deg = lat_AA_deg, time_of_day = time_of_days[i], days_after_vernal_equinox = days_after_vernal_equinoxs[i], longitude_deg = long_AA_deg, UTC_offset = 10))
    print("")


print("\n\n\n")

print(horizontal_from_equatorial(declination = rad_from_degree_min_sec("48:19:27.7"), RA = radians_from_HHMMSS("19:15:23.63"), latitude_deg = lat_McDonald_deg, time_of_day = "00:23:49.4", days_after_vernal_equinox = 90, longitude_deg = long_McDonald_deg, UTC_offset = -5))



print("\n\n")


LST_McDonald = get_LST(local_time_of_day = "11:23:19.2", days_after_vernal_equinox = 89, longitude_deg = long_McDonald_deg, UTC_offset = -5)

print("LST at McDonald Obs: "+str(LST_McDonald))



LST_AA = get_LST(local_time_of_day = "02:23:19.2", days_after_vernal_equinox = 90, longitude_deg = long_AA_deg, UTC_offset = 10)

print("LST at AA Obs: "+str(LST_AA))




print("LST proper: "+str(get_LST_good(local_time = "13:41:25", JD = 2459747.5, longitude_deg = -77.8600012, UTC_offset = -4)))
