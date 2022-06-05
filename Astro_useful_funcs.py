import numpy as np
import math

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
       
    if num_to_clean > 1:
        #looks like 457832948598 or maybe like 54783928457394.0000004 or like 548372945833.45387438532
        #either way, it needs to look like 5.6584E16 or whatever-specifically, with only clean_arg-1 digits after the decimal
        
        power = math.floor(np.log10(num_to_clean))
        C = round(num_to_clean/10.**power, num_digits-1)
        cleaned_num = float(str(C)+"E"+str(power))
    else:
        #looks like 0.0438298 or 0.33333333331 or whatever. This bit is easy, we just need to round the decimal
        cleaned_num = round(num_to_clean, num_digits-1)
    return cleaned_num
    
def str_clean_number(num_to_clean, num_digits = 5, max_power = 4):
    '''inputs: num_to_clean-the number to be cleaned and be made to look pretty, num_digits-default to 5, if passed can set how many digits to have in the number, max_power-default to 3, if passed can set maximum power to display without scientific notation
       outputs: cleaned_num-number that has been processed and will be returned
       this function takes in a number that may be very large (13678323450) or very small (0.00000000554) or have some extra small digit added to it (4534.00000001) and will return a pretty version of the same number as a string, so it can be displayed'''
       
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
