import numpy as np
import math


G_SI = 6.674E-11 #Universal gravitational constant G in SI units N*m^2/kg^2
G_cgs = 6.674E-8 #Universal gravitational constant G in cgs units dyne*cm^2/g^2

H_0_km_s_Mpc = 67.37 #hubble constant in units of km/s/Mpc
H_0_hz = 2.1833E-18 #hubble constant in units of 1/sec, this is useful for SI unit calculations
H_0 = H_0_km_s_Mpc
H_0_SI = H_0_hz

h = H_0_km_s_Mpc/100 #dimensionless h constant for Hubble constant
rho_crit_SI = H_0_SI**2*3/(8*np.pi*G_SI) #critical density of the universe in SI units
rho_crit_cgs = 1.9E-29*h #in cgs units of g/cm3

Omega_M = 0.3147 #mass density of the universe
Omega_Lambda = 0.6847 #dark energy density of universe
Omega_k = 0.001
t_0_Gyr = 13.797 #age of universe in Gyr
t_0 = t_0_Gyr


def find_nearest(wv,value):
    '''inputs: wv-list of data to look in, value-the value we are looking for in wv
       outputs: idx-index in wv that is closest to value
       this finds the index of the closest element to value in list wv'''
    idx= (np.abs(np.asarray(wv) - value)).argmin()
    return idx
    
def gaussian(x, mu, sig, A = 1.0):
    '''inputs: x-list or array of independent variables over which to build a gaussian, mu-center of gaussian, sig-standard deviation of gaussian, A-vertical scaling constant of gaussian (default to 1.0)
       outputs: gaussian-array of dependent variable defining gaussian with input parameters
       This function passes out a gaussian with the parameters that are passed in'''
    return np.asarray(A*np.exp(-np.power(np.asarray(x) - mu, 2.) / (2 * np.power(sig, 2.))))
    
def correct_cosmological_redshift(wavelengths, fluxes, z):
    '''inputs: wavelengths-list of wavelengths for the spectrum that is passed in, fluxes-corresponding list of fluxes for the spectrum, z-cosmological redshift of the obeject emitting this spectrum
       outputs: wavelengths_prog_frame-list of wavelengths of the spectrum in progenitor frame, fluxes_prog_frame-list of corresponding fluxes of the spectrum in progenitor frame
       This function takes an observed spectrum and converts it into a progenitor frame spectrum, given a redshift z'''
    fluxes_prog_frame = [i*(1+z) for i in fluxes]
    wavelengths_prog_frame = [i/(1+z) for i in wavelengths]
    
    return wavelengths_prog_frame, fluxes_prog_frame

def power_from_obs_flux(obs_flux, distance):
    '''inputs: obs_flux-the observed flux observation of an object, distance-distance to the object
       outputs: power-total power output of object
       This function gives the total power output of an object in units of [obs_flux]*[distance]^2--as in, whatever units are passed in are the units that come out. Note that this just gives the total power output, not the luminosity at some distance or some measure of apparent magnitude or absolute magnitude. Also note that this function inherently assumes this object is an isotropic emitter.'''
    
    power = obs_flux*4*np.pi*distance**2.0
    return power
    


def get_absolute_mag(app_mag, dist_parsec):
	'''inputs: app_mag-apparent magnitude of an object, dist_parsec-distance to the object in parsecs
	   outputs: absolute_mag-absolute magnitude of the object, dist_mod-distance modulus of the object
	   this function uses the equation m-M = 5log(d/10 parsec) to get the distance modulus and the absolute magnitude of an object, given its apparent magnitude and distance'''
	
	dist_mod = 5*np.log10(dist_parsec/10)
	absolute_mag = -1.0*dist_mod - app_mag
	
	return {"absolute_mag": absolute_mag, "dist_mod":dist_mod}


def fnu_from_flambda(flambda, lam, c = 3.0E18):
	'''inputs: flambda-f_lambda that we want to get the f_nu of, lam-wavelength at which this is happening (angstroms is default), c-optional speed of light (this is useful to specify units of length, if something other than angstroms is needed, specify speed of light in these length units)
	   outputs: fnu-f_nu from converting flambda above, units are returned in same units of flambda as passed in, multiplied by length units passed in divided by seconds
	   this converts flambda into fnu for a given measurement. The default length units are angstroms and can be changed by changing the units of c, default times units are seconds, and cannot be changed. The function returns fnu in the units of flambda multiplied by the length units specified and divided 
by seconds'''

	fnu = flambda*lam**2/c
	return fnu


def flambda_from_fnu(fnu, lam, c = 3.0E18):
    '''inputs: fnu-f_nu that we want to get the f_lambda of, lam-wavelength at which this is happening (angstroms is default), c-optional speed of light (this is useful to specify units of length, if something other than angstroms is needed, specify speed of light in these length units)
       outputs: flambda-f_lambda from converting nuabove, units are returned in same units of flambda as passed in, multiplied by seconds divided by length units passed in
       this converts flambda into fnu for a given measurement. The default length units are angstroms and can be changed by changing the units of c, default times units are seconds, and cannot be changed. The function returns flambda in the units of fnu multiplied by seconds divided by the length units passed in'''

    flambda = fnu*c/lam**2
    return flambda


def flux_nu_jansky_from_flux_nu_cgs(flux_nu_cgs):
    '''inputs: flux_nu_cgs-flux in units of erg/cm^2/s/Hz
       outputs: flux_nu_jansky-flux in units of janskies
    '''
    return flux_nu_cgs * 10.**23


def get_AB_mag_from_flux_nu(flux_nu_cgs = -np.pi, flux_nu_jansky = -np.pi):
    '''inputs: flux_nu_cgs-flux in units of erg/cm^2/s/Hz, flux_nu_jansky-flux in units of janskies
       outputs: AB_mag-AB magnitude measure of magnitude from flux that is passed in
       This function gets the AB magnitude of a f_nu measurement. Either flux_nu_cgs or flux_nu_jansky needs to be passed in (or both). If neither is passed in, an error is raised
    '''
    if flux_nu_cgs == -np.pi and flux_nu_jansky == -np.pi:
        print("Error: flux not passed in")
        return np.pi
    elif flux_nu_cgs != -np.pi and flux_nu_jansky != -np.pi and flux_nu_jansky_from_flux_nu_cgs(flux_nu_cgs) != flux_nu_jansky:
        print("Error: fluxes not consistent")
        return np.pi
    elif flux_nu_cgs != -np.pi:
        flux_nu_jansky = flux_nu_jansky_from_flux_nu_cgs(flux_nu_cgs)
    
    AB_mag = -2.5*np.log10(flux_nu_jansky)+8.9
    
    return AB_mag
    
def get_AB_mag_from_flux_lam(lam, flux_lam_cgs = -np.pi):
    '''inputs: lam-wavelength of this f_lam measurement in Angstroms, flux_lam_cgs-flux in units of erg/cm^2/s/angstrom
       outputs: AB_mag-AB magnitude measure of magnitude from flux that is passed in
       This function gets the AB magnitude of a f_lam measurement using the lam that is passed in. Note that AB_mag is defined natively in f_nu space, so a wavelength must be given to properly convert to f_nu space. This turns f_lam into f_nu_jansky, then gets AB-mag from that
    '''
    if flux_lam_cgs == -np.pi:
        print("Error: flux not passed in")
        return np.pi

    flux_nu_cgs = fnu_from_flambda(flambda = flux_lam_cgs, lam = lam)
    flux_nu_jansky = flux_nu_jansky_from_flux_nu_cgs(flux_nu_cgs = flux_nu_cgs)
    AB_mag = get_AB_mag_from_flux_nu(flux_nu_jansky = flux_nu_jansky)
    
    return AB_mag
        
def get_tangential_vel(proper_motion, dist_parsec):
    '''inputs: proper_motion-proper motion of object in arsec/yr, dist_parsec-distance to object in parsec
       outputs: tang_vel-tangential velocity of object in km/sec
       this gets the tangential velocity of an object in km/sec, given its proper motion in arcsec/yr and distance in parsecs
    '''
    tang_vel = 4.74*proper_motion*dist_parsec
    return tang_vel
   
def get_dist_tang_vel_method(radial_vel, lam_rad, proper_motion):
    '''inputs: radial_vel-radial velocity of object in km/sec, lam_rad between direction of radial velocity and direction of overal motion in radians, proper_motion-proper motion of object in arcsec/yr
       outputs: dist-distance to object in parsec
       this gets the distance to an object using the moving cluster method, given the radial velocity, angle between radial velocity and overall convergence direction, and proper motion
    '''
    
    dist = radial_vel*np.tan(lam_rad)/(4.74*proper_motion)
    
    return dist

def radian_from_degree(angle_degree):
    '''inputs: angle_degree-an angle measured in degrees
       outputs: angle_rad-angle measured in radians
       gets radians from degrees of an angle
    '''
    return np.pi*angle_degree/180
    
def deg_min_sec_from_degree(degree):
    '''inputs: degree-angle measured in degrees
       outputs: deg_min_sec-angle measured in units of degrees, minutes, seconds
       This function converts an angle in units of degrees to degrees minutes and seconds'''
       
    deg = int(degree)
    
    deg_frac = degree-deg
    minutes = deg_frac*60.
    
    minutes_frac = minutes-int(minutes)
    secs = minutes_frac*60.
    
    if deg < 0:
        minutes *= -1.0
        secs *= -1.0
    
    return str(deg)+":"+str(int(minutes))+":"+str(round(secs, 4))
    
    
    
    
def rad_from_degree_min_sec(angle_degree_min_sec):
    '''inputs: angle_degree_min_sec-angle measured in units of degrees, arcminutes, and arcseconds
       outputs: angle_rad-corresponding angle in radians
       This function gets radians from an angle given in units of degreesArcminArcsec'''
       
    GG = angle_degree_min_sec.split(":")
    deg_floor = float(GG[0])
    arc_min = float(GG[1])
    arc_sec = float(GG[2])
    
    if deg_floor < 0:
        arc_min *= -1.0
        arc_sec *= -1.0

    angle_deg = deg_floor + arc_min/60.+arc_sec/3600.
    
    return radian_from_degree(angle_deg)

def energy_from_wavelength(lam_m):
    '''inputs: lam_m-wavelength in METERS that we want the energy of
       outputs: lam_nm-wavelength of light in nanometers, lam_m-wavelength of light in meters, lam_a-wavelength of light in angstroms, energy_erg-energy of this light in ergs, energy_ev-energy of this light in ev, energy_joule-energy of this light in joules, freq-frequency of this light in hz, correct_color-string representing what band this would be (like 'radio', 'visible', 'gamma' or whatever), t_eff-effective temperature of this light, assume h*frequency = boltzmann constant*t_eff
       This function takes a wavelength in meters and gives back all the info of this light-including the energy in various units, frequency, effective temperature, and where on the EM spectrum it lies
    '''
    c_m = 3.0E8 #speed of light in meters/sec
    freq = c_m/lam_m #frequency of light in hz
    
    h = 6.626E-34 #planck constant in SI units
    energy_joule = h*freq #energy of light in Joules
    
    energy_ev = energy_joule*6.242E18 #energy of light in ev
    energy_erg = energy_joule*1.0E7 #energy of light in ergs
    
    kB_SI = 1.381E-23 #boltzmann constant in SI units
    t_eff = energy_joule/kB_SI #effective temperature of this light in Kelvin
    
    spectrum_wavelength_bounds = {"Gamma": [0, 1.0E-12], "Hard X-Ray": [1.0E-12, 1.0E-10], "Soft X-Ray": [1.0E-10, 3E-9], "Extreme Ultraviolet": [3E-9, 30E-9], "Near Ultraviolet": [30E-9, 380E-9], "Violet": [380E-9, 450E-9], "Blue":[450E-9, 485E-9], "Cyan":[485E-9, 500E-9], "Green":[500E-9, 565E-9], "Yellow":[565E-9, 590E-9], "Orange":[590E-9, 625E-9], "Red":[625E-9, 750E-9], "Near Infrared":[750E-9, 3E-6], "Mid Infrared":[3E-6,15E-6], "Far Infrared":[15E-6, 1E-3], "Microwave (EHF)":[1E-3, 10E-3], "Microwave (SHF)":[1E-2, 10E-2], "Microwave (UHF)":[1E-1, 1], "Radio (VHF)":[1, 10], "Radio (HF)":[10, 100], "Radio (MF)":[100, 1000], "Radio (LF)":[1E3, 10E3], "Radio (VLF)":[10E3, 100E3], "Radio (ULF)":[100E3, 1000E3], "Radio (SLF)":[1.0E6, 1.0E7], "Radio (ELF)":[1.0E7, 1.0E8]}
    
    correct_color = 'Not found'
    for color in spectrum_wavelength_bounds.keys():
        color_bounds = spectrum_wavelength_bounds[color]
        
        if color_bounds[0] <= lam_m and lam_m < color_bounds[1]:
            correct_color = color
    

    return {"lam_nm":lam_m*1.0E9, "lam_m":lam_m, "lam_a":lam_m*1.0E10, "energy_erg":energy_erg, "energy_ev":energy_ev, "energy_joule":energy_joule, "freq": freq, "correct_color":correct_color, "t_eff":t_eff}

def wavelength_from_energy(energy_joule):
    '''inputs: energy_joule-energy of light in joules
       outputs: lam_nm-wavelength of light in nanometers, lam_m-wavelength of light in meters, lam_a-wavelength of light in angstroms, energy_erg-energy of this light in ergs, energy_ev-energy of this light in ev, energy_joule-energy of this light in joules, freq-frequency of this light in hz, correct_color-string representing what band this would be (like 'radio', 'visible', 'gamma' or whatever), t_eff-effective temperature of this light, assume h*frequency = boltzmann constant*t_eff
       This function takes an energy in joules and gives back all the info of this light-including the energy in various units, frequency, effective temperature, and where on the EM spectrum it lies
    '''
    h = 6.626E-34 #planck constant in SI units
    c = 3.0E8 #speed of light in SI units
    lam_m = c/(energy_joule/h)
    
    d = energy_from_wavelength(lam_m = lam_m)
    
    return d
    
def wavelength_from_temp(temp_k): #this function will be overloaded with different measures of energy, to turn into wavelength
    '''inputs: temp_k-temperature in Kelvin
       outputs: lam_nm-wavelength of light in nanometers, lam_m-wavelength of light in meters, lam_a-wavelength of light in angstroms, energy_erg-energy of this light in ergs, energy_ev-energy of this light in ev, energy_joule-energy of this light in joules, freq-frequency of this light in hz, correct_color-string representing what band this would be (like 'radio', 'visible', 'gamma' or whatever), t_eff-effective temperature of this light, assume h*frequency = boltzmann constant*t_eff
       This function takes a temperature in Kelvin and gives back all the info of this light-including the energy in various units, frequency, effective temperature, and where on the EM spectrum it lies
    '''
    kB_SI = 1.381E-23 #boltzmann constant in SI units
    d = wavelength_from_energy(energy_joule = kB_SI*temp_k)
    
    
    return d
    
    
def clean_number(num_to_clean, num_digits = 5):
    '''inputs: num_to_clean-the number to be cleaned and be made to look pretty, clean_arg-default to 5, if passed can set how many digits to have in the number
       outputs: cleaned_num-number that has been processed and will be returned
       this function takes in a number that may be very large (13678323450) or very small (0.00000000554) or have some extra small digit added to it (4534.00000001) and will return a pretty version of the same number as a float'''
    try:
        num_to_clean = float(num_to_clean)
    except ValueError as e:
        return num_to_clean
    
    if num_to_clean > 1:
        #looks like 457832948598 or maybe like 54783928457394.0000004 or like 548372945833.45387438532
        #either way, it needs to look like 5.6584E16 or whatever-specifically, with only clean_arg-1 digits after the decimal
        
        power = math.floor(np.log10(num_to_clean))
        C = round(num_to_clean/10.**power, num_digits-1)
        cleaned_num = float(str(C)+"E"+str(power))
    else:
        #looks like 0.0438298 or 0.33333333331 or whatever. This bit is easy, we just need to round the decimal, if the rounded version is just zero, we'll have to mess around a little
        cleaned_num = round(num_to_clean, num_digits-1)
    return cleaned_num
    
def str_clean_number(num_to_clean, num_digits = 5, max_power = 4):
    '''inputs: num_to_clean-the number to be cleaned and be made to look pretty, num_digits-default to 5, if passed can set how many digits to have in the number, max_power-default to 3, if passed can set maximum power to display without scientific notation
       outputs: cleaned_num-number that has been processed and will be returned
       this function takes in a number that may be very large (13678323450) or very small (0.00000000554) or have some extra small digit added to it (4534.00000001) and will return a pretty version of the same number as a string, so it can be displayed'''
    try:
        num_to_clean = float(num_to_clean)
    except ValueError as e:
        return num_to_clean
    
    
    if num_to_clean == 0.0:
        return "0.0"
    if num_to_clean > 1:
        #looks like 457832948598 or maybe like 54783928457394.0000004 or like 548372945833.45387438532
        #either way, it needs to look like 5.6584E16 or whatever-specifically, with only clean_arg-1 digits after the decimal
        
        power = math.floor(np.log10(num_to_clean))
        C = round(num_to_clean/10.**power, num_digits-1)
        if power <= max_power:
            cleaned_num = str(C*10.**power)
        else:
            cleaned_num = str(C)+"E"+str(power)
    else:
        #looks like 0.0438298 or 0.33333333331 or whatever. This bit is easy, we just need to round the decimal
        cleaned_num = str(round(num_to_clean, num_digits-1))
        
        
        if cleaned_num == "0.0" and num_to_clean != 0.0:
            cleaned_num = str(num_to_clean)[0:num_digits]
            power = math.floor(np.log10(num_to_clean))
            
            C = round(num_to_clean/10.**power, num_digits-1)
            
            cleaned_num = str(C)+"E"+str(power)
    return cleaned_num
    
#reddening, de-reddening
def get_reddening_extinction(R_lam, E_B_V):
    
    '''inputs: R_lam-the R value for reddening in A_lam = R_lam*E(B-V), E_B_V-the reddening value E(B-V)
       outputs: A_lam-the extinction in magnitudes using A_lam = R_lam*E(B-V) at a given wavelength lam
       This function just gives the number of magnitudes that are extinguished by reddening'''
    A_lam = R_lam*E_B_V
    
    return A_lam
    
def get_E_B_V(H_alpha_H_beta_flux_obs_ratio, H_alpha_H_beta_flux_true_ratio = 2.86, R_HBeta = 3.609, R_HAlpha = 2.535):
    '''inputs: H_alpha_H_beta_flux_obs_ratio-ratio of the flux from H-alpha to the flux from H-beta, H_alpha_H_beta_flux_true_ratio-the true value of this flux, this is a constant known to be 2.86, R_HBeta-the R value for H_Beta (default to 3.609, the avg for Milky Way), R_HAlpha-the R value for H_Alpha (default to 2.535, the avg for Milky Way)
       outputs: E_B_V-value of E(B-V) for the given flux ratio
       this function assumes Milky Way values for reddening of both H-Alpha and H-Beta and gives the corresponding E_B_V'''
       
    E_B_V = 2.5*(np.log10(H_alpha_H_beta_flux_obs_ratio)-np.log10(H_alpha_H_beta_flux_true_ratio))/(R_HBeta-R_HAlpha)
    
    return E_B_V

def deredden_flux(f_lam_obs, R_lam, H_alpha_H_beta_flux_obs_ratio, H_alpha_H_beta_flux_true_ratio = 2.86, R_HBeta = 3.609, R_HAlpha = 2.535):
    '''inputs: f_lam_obs-observed (reddened) flux at a given wavelength, H_alpha_H_beta_flux_obs_ratio-ratio of the flux from H-alpha to the flux from H-beta, H_alpha_H_beta_flux_true_ratio-the true value of this flux, this is a constant known to be 2.86, R_HBeta-the R value for H_Beta (default to 3.609, the avg for Milky Way), R_HAlpha-the R value for H_Alpha (default to 2.535, the avg for Milky Way)
       outputs: f_lam_true-true (dereddened) value of f_lam
       this function assumes Milky Way values for reddening of both H-Alpha and H-Beta and finds the de-reddened flux for a given observed flux'''
    
    E_B_V = get_E_B_V(H_alpha_H_beta_flux_obs_ratio = H_alpha_H_beta_flux_obs_ratio, H_alpha_H_beta_flux_true_ratio = 2.86, R_HBeta = 3.609, R_HAlpha = 2.535)
    
    
    f_lam_true = f_lam_obs*10.**(R_lam*E_B_V/2.5)
    return f_lam_true
    
       
#density diagnostics

def density_diagnostic(f_3729, f_3726):
    '''inputs: f_3729-flux in 3729 angstroms, f_3726-flux in 3726 angstroms
       outputs: density-gives density in units of Ne/sqrt(Te/10^4), judgement-string that says it's low, medium, or high
       This function uses the figure in slide deck 15b from Robin's lecture series given the two Oxygen line fluxes to get whether this density is low, medium, or high, if it's low, we can use temperature diagnostics'''
    
    density_data = np.loadtxt("Figures/density_diagnostic_3729_3726_ratio.txt")
    dens = density_data[:,0]
    intensity = data[:,1]
    idx = find_nearest(wv = intensity, value = f_3729/f_3726)
    
    density = dens[idx]
    
    judgement = ""
    if dens < 100:
        judgement = "low"
    elif dens < 3.0E3:
        judgement = "med"
    else:
        judgement = "high"
    
    return {"density": density, "judgement": judgement}
    
#temperature diagnostics
def temp_diagnostic_low_density(f_5007, f_4959, f_4363):
    '''inputs: f_5007-flux in 5007 angstroms, f_4959-flux in 4959 angstroms, f_4363-flux in 4363 angstroms
       outputs: temp-temperature found from plot
       This function uses the figure in slide deck 15b from Robin's lecture series to get temperature given these three Oxygen line fluxes. Remember that this only works if the Ne density is low-check that first!'''
       
    temp_data = np.loadtxt("Figures/temp_diagnostic_5007_4959_4363.txt")
    temps = temp_data[:,0]
    ratios = temp_data[:,1]
    idx = find_nearest(wv=ratios, value = (f_5007+f_4959)/f_4363)
    
    temp = temps[idx]
    
    return temp


def get_temp_density_diagnostic(f_5007, f_4959, f_4363, f_3729, f_3726):
    '''inputs: f_3729-flux in 3729 angstroms, f_3726-flux in 3726 angstroms, f_5007-flux in 5007 angstroms, f_4959-flux in 4959 angstroms, f_4363-flux in 4363 angstroms
       outputs: temp-temperature found, density-Ne in cm^-3
       This function uses the figures in slide deck 15b from Robin's lecture series to get density using the two oxygen forbidden lines 3729 and 3726, and if density is low enough, then uses the Oxygen forbidden lines 5007, 4959, and 4363 to get temperature'''
    dd = density_diagnostic(f_3729 = f_3729, f_3726 = f_3726)
    density_odd_units = dd['density']
    judgement = dd['judgement']
    
    if judgement == 'low':
        temp_K = temp_diagnostic_low_density(f_5007 = f_5007, f_4959 = f_4959, f_4363 = f_4363)
        density_cm3 = density_odd_units*np.sqrt(temp_K/10000)
        return {"temp_K": temp_K, "density_cm3": density_cm3}
    
    return {"density_odd_units": density_odd_units}

#cosmological distance measurements

def physical_dist_from_comoving_dist(comoving_dist, scale_factor):
    '''inputs: comoving_dist-distance to object in comoving coordinates, scale_factor-scale factor of universe
       outputs: physical_dist-physical distance to object (units are same as units in comoving_dist
       this function takes in a comoving distance and gives a physical distance to an object. The units of the answer are the same as would come from comoving_dist'''
       
    
    return comoving_dist*scale_factor
    
def z_from_scale_factor(a, a_0=1.0):
    '''inputs: a-scale factor of universe at some time t, a_0-(optional) scale factor of universe at time 0-defaults to t = now and a_0=1.0
       outputs: z-redshift measurement of the object at time t, when seen from time 0
       This function gets the redhisft measurement of an object given the scale factor of the universe at that time. This function is flexible enough to get the redshift that one object in the past would see another object in the further past. We use the equation (1+z) = a_0/a'''
    
    return 1-(a_0/a)
    
def scale_factor_from_z(z, a_0 = 1.0):
    '''inputs: z-redshift of an object in the distance, a_0-(optional) scale factor of the universe when the object was observed-defaults to a_0=1.0 corresponding to now
       outputs: a-scale factor of the universe when this object emitted the light being observed
       This function turns a redshift into a scale factor using the equation (1+z) = a_0/a'''

    return a_0/(1+z)
    
#lookback time as a function of z, as a function of a

def conv_time_sec(time_sec):
    '''inputs: time_sec-time measured in seconds
       outputs: time_dic-dictionary of time_sec converted into various units, including minutes, hours, solar days, years, Myr, and Gyr
       This function converts time in seconds to a bunch of other units, returning all of them in a dictionary'''
    
    time_min = time_sec/60.0
    time_hr = time_sec/3600.0
    time_day = time_sec/86400.0
    time_yr = time_day/365.24219
    time_Myr = time_yr/1.0E6
    time_Gyr = time_yr/1.0E9
    
    return {"time_sec":time_sec,"time_min":time_min, "time_hr":time_hr, "time_day":time_day, "time_yr": time_yr, "time_Myr":time_Myr, "time_Gyr":time_Gyr}



def milne_lookback_time(z, hubble_const = H_0_SI):
    '''inputs: z-redshift of some distant object, hubble_const-(optional) current hubble constant of the universe-defaults to the value defined above
       outputs: delta_t-dictionary of lookback time measurment from now (t_0) to some object at redshift z in various units
       This function gets the lookback time to an object at redshift z using the equation t = (1/H_0)(1-1/(1+z)). Note that this assumes the Milne cosmology, which says that \rho = 0 and that the universe is empty of matter'''
    
    lookback_time_sec = (1.0/hubble_const)*(1.0-1.0/(1.0+z))
    
    time_dic = conv_time_sec(time_sec = lookback_time_sec)
    lookback_time_yr = time_dic['time_yr']
    lookback_time_Gyr = time_dic['time_Gyr']
    return {"lookback_time_sec":lookback_time_sec, "lookback_time_yr":lookback_time_yr, "lookback_time_Gyr":lookback_time_Gyr}



def de_sitter_lookback_time(z, hubble_const = H_0_SI):
    '''inputs: z-redshift of some distant object, hubble_const-(optional) current hubble constant of the universe-defaults to the value defined above
       outputs: delta_t-dictionary of lookback time measurment from now (t_0) to some object at redshift z in various units
       This function gets the lookback time to an object at redshift z using the equation t = (2/3)(1/H_0)(1-(1+z)^(-3/2)). Note that this assumes the Einstein-de Sitter cosmology, which says that E = 0 and the universe is a critical universe, where total energy is zero'''
       
    lookback_time_sec = (2.0/3.0)*(1.0/hubble_const)*(1.0-(1.0+z)**(-1.5))
    
    time_dic = conv_time_sec(time_sec = lookback_time_sec)
    lookback_time_yr = time_dic['time_yr']
    lookback_time_Gyr = time_dic['time_Gyr']
    return {"lookback_time_sec":lookback_time_sec, "lookback_time_yr":lookback_time_yr, "lookback_time_Gyr":lookback_time_Gyr}

def radians_from_HHMMSS(angle_HHMMSS):
    '''inputs: angle_HHMMSS-angle in units of HH:MM:SS
       outputs: angle_rad-the angle passed in units of radians
       This function turns an angle in units of hours, minutes, seconds into an angle in radians'''
       
    g = angle_HHMMSS.split(":")
    HH = float(g[0])
    
    
    MM = float(g[1])
    SS = float(g[2])
    
    if HH < 0:
        MM *= -1.0
        SS *= -1.0
    
    #print("HH: "+str(HH))
    #print("MM: "+str(MM))
    #print("SS: "+str(SS))
    
    angle_frac = HH/24. + MM/(24.*60.) + SS/(24.*60.*60.)
    angle_rad = angle_frac * 2.*np.pi
    
    return angle_rad
       
    


def HHMMSS_from_radians(angle_rad):
    '''inputs: angle_rad-the angle passed in units of radians
       outputs: angle_HHMMSS-angle in units of HH:MM:SS
       This function turns an angle in radians into an angle in units of hours, minutes, seconds'''


    if angle_rad < 0:
        angle_rad += 2*np.pi
        
        
    angle_frac = angle_rad/(2.*np.pi)
    
    hours = angle_frac*24.
    frac_hours = hours - math.floor(hours)
    
    mins = frac_hours * 60.
    frac_mins = mins - math.floor(mins)
    
    secs = frac_mins * 60.
    
    return str(math.floor(hours))+":"+str(math.floor(mins))+":"+str(round(secs, 5))


def equatorial_from_horizontal(altitude, azimuth, latitude_deg, time_of_day, days_after_vernal_equinox, longitude_deg, UTC_offset):
    '''inputs: altitude-altitude angle of some object in radians, azimuth-azimuth angle of some object in radians, latitude-latitude of location on Earth in radians
       outputs: RA-right ascension of the object in HH:MM:SS, dec-declination of the object in radians
       This function takes altitude, azimuth, and latitude to give the right ascension and declination of some object measured on the sky plane.'''
       
    
    latitude = latitude_deg * np.pi/180.
    #print("altitude: "+str(altitude))
    #print("azimuth: "+str(azimuth))
    #print("lat: "+str(latitude))
    #print("time of day: "+str(time_of_day))
    #print("days_after..: "+str(days_after_vernal_equinox))
    
    current_LST_hrs = get_LST(local_time_of_day = time_of_day, days_after_vernal_equinox = days_after_vernal_equinox, longitude_deg=longitude_deg, UTC_offset = UTC_offset)['LST_hrs']
    
    
    #print("current_LST_hrs: "+str(current_LST_hrs))
    
    current_LST_radians = current_LST_hrs/24.*2*np.pi
    
    #print("current_LST_radians: "+str(current_LST_radians))
    
    #we'll use the equation sin(dec) = -cos(azimuth)*cos(altitude)*cos(lat) + sin(altitude)*sin(lat) to get declination first
    
    dec = math.asin(math.cos(azimuth)*math.cos(altitude)*math.cos(latitude)+math.sin(altitude)*math.sin(latitude)) #declination can only be between -np.pi/2 to np.pi/2, so there's no degeneracy in this function
    
    #print("dec: "+str(dec))
    #print("")
    

    #I need to use another equation to tell which of h1 or h2 is the correct right ascension
    #I can use sin(h)*cos(dec) = -sin(azimuth)*cos(altitude) to get h3, the asin(..) solution
    
    h3 = math.asin(-1.0*math.sin(azimuth)*math.cos(altitude)/math.cos(dec))

    
    h_radians = h3
    #print("found h_radians: "+str(h_radians)+" but I wanted to use: "+str(A[combo[0]]))
    
    RA_radians = current_LST_radians-h_radians
    RA_HHMMSS = HHMMSS_from_radians(RA_radians)
    
    dec_rad = dec
    dec_HHMMSS = HHMMSS_from_radians(dec_rad)
    dec_deg = dec_rad*180./np.pi
    dec_deg_m_s = deg_min_sec_from_degree(dec_deg)
    
    RA = RA_HHMMSS
    dec = dec_deg_m_s
    
    return {"RA": RA, "dec": dec}
    
    
    
    
def get_UTC_offset(longitude_deg):
    '''inputs: longitude_deg-longitude of location in degrees (negative for west, positive for east)
       outputs: UTC_offset-an integer telling how many hours to add to local time at this longitude to get UTC time
       This function gives how many hours behind the longitude is from UTC time'''
    from datetime import datetime
    datetime.utcnow()
    
    hours_offset = math.round(longitude/15)
    
    return hours_offset #negative 1 because we want how many hours <<to add>>
    

def get_LST(local_time_of_day, days_after_vernal_equinox, longitude_deg, UTC_offset = -45):
    '''inputs: time_of_day-the time of day in format HH:MM:SS as a string, days_after_vernal_equinox-days since vernal equinox for the date in question, longitude_deg-the longitude of the location on Earth measured in degrees (negative for West and positive for East), UTC_offset-(optional) the number of hours to add to UTC to get local time this defaults to -45 and if it is set to -45 this will attempt to calculate UTC time using the longitude
       outputs: LST-corresponding local sidereal time of this time measure
       This function gets local sidereal time for a given time measurement on a given day. The equation I am using is LST = T + 12hours + (days after vernal equinox)*4 min'''
    
    GG_TOD = local_time_of_day.split(":")
    local_time = float(GG_TOD[0])+float(GG_TOD[1])/60.+float(GG_TOD[2])/3600.
    
    #print("local time: "+str(local_time))
    #1) go from local time to UTC
    
    if UTC_offset == -45:
        UTC_offset = get_UTC_offset(longitude_deg = longitude_deg)
    
    UTC_time = local_time-UTC_offset
    
    #print("UTC time: "+str(UTC_time))
    
    
    #2) get GST-sidereal time at the Greenwich meridian
    
    GST_hrs = (UTC_time + 12 + days_after_vernal_equinox*(3./60.+55./3600.))%24
    #print("GST: "+str(GST_hrs))
    
    
    #3) convert from GST to LST using the local longitude in HHMMSS format
    
    longitude_hrs = longitude_deg/360. * 24.
    LST_hrs = GST_hrs+longitude_hrs
    
    
    
    #LST_hrs = time_of_day_hrs + 12 + days_after_vernal_equinox*4./60.
    #LST_hrs = (time_of_day_hrs + 12 + days_after_vernal_equinox*3.9333/60.)%24.
    
    LST_hr_floor = math.floor(LST_hrs)
    
    LST_mins = (LST_hrs - LST_hr_floor)*60.
    LST_mins_floor = math.floor(LST_mins)
    
    LST_secs = round((LST_mins - LST_mins_floor)*60.,4)
    
    LST_HHMMSS = str(LST_hr_floor)+":"+str(LST_mins_floor)+":"+str(LST_secs)
    
    return {"LST_HHMMSS": LST_HHMMSS, "LST_hrs": LST_hrs}
    
    
def get_LST_good(local_time, JD, longitude_deg, UTC_offset):
    
    GG_TOD = local_time.split(":")
    local_time = float(GG_TOD[0])+float(GG_TOD[1])/60.+float(GG_TOD[2])/3600.
    
    #print("local time: "+str(local_time))
    #1) go from local time to UTC
    
    if UTC_offset == -45:
        UTC_offset = get_UTC_offset(longitude_deg = longitude_deg)
    
    UTC_time = local_time-UTC_offset
    
    T = (JD - 2451545)/36525
    theta = 280.46061837 + 360.98564736629 * (JD - 2451545.0) + (0.000387933 * T * T) - (T * T * T / 38710000.0)
    
    A = theta*360/(2*np.pi) + longitude_deg
    print("A: "+str(A*24/360 % 24))
    
    return theta%24
    
    

def horizontal_from_equatorial(declination, RA, latitude_deg, time_of_day, days_after_vernal_equinox, longitude_deg, UTC_offset):
    #a is altitude, A is azimuth
    
    latitude = latitude_deg * np.pi/180.
    
    current_LST_hrs = get_LST(local_time_of_day = time_of_day, days_after_vernal_equinox = days_after_vernal_equinox, longitude_deg=longitude_deg, UTC_offset = UTC_offset)['LST_hrs']
    
    current_LST_radians = (current_LST_hrs/24.)*2.*np.pi
    #print("declination: "+str(declination))
    #print("RA: "+str(RA))
    #print("current_LST_hrs: "+str(current_LST_hrs))
    
    h = current_LST_radians-RA #this works a good bit better if I subtract 0.5 from h
    
    
    
    
    #print("current_LST_radians:" +str(current_LST_radians))
    #print("RA: " +str(RA))
    #print("h: "+str(h))
    
    
    a = math.asin(math.cos(h)*math.cos(declination)*math.cos(latitude)+math.sin(declination)*math.sin(latitude))
    
    #use cos(a)*sin(A) = -cos(delta)*sin(H)
    A = math.asin(-1.0*math.cos(declination)*math.sin(h)/math.cos(a))
    
    if A < 0:
        A += 2.*np.pi
    
    return{"alt": deg_min_sec_from_degree(a*180./np.pi), "az":deg_min_sec_from_degree(A*180./np.pi)}

