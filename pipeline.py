import os
import shutil
from astropy.io import fits
import sys
import numpy as np
import logging
import multiprocessing
from multiprocessing import Pool
import csv
import pandas as pd
import matplotlib.pyplot as plt
import math 
from statistics import mean
import ipywidgets as widgets

from ipywidgets import HBox, VBox


from IPython.display import display


from setup import datapathWidget,number,coadd, res_widget,ispecpathWidget, input_widget, osx

ispec_dir = ispecpathWidget.value
if ispec_dir[-1] != '/':
            ispec_dir=ispec_dir+'/'
import ispec 
oss = osx.value
if oss == 'Mac/iOS':
    code = 'moog'
else:
    code = 'width'
    
    
def parameters_from_ew(b):
    
    
    
    resolution=res_widget.value
    mask = input_widget.value

    to_resolution = resolution
    if number.value=='Multiple':
        paths = [ f.path for f in os.scandir(datapathWidget.value) if f.is_dir()]
    if number.value=='Single':
        paths = [datapathWidget.value]            

    file_name = '/prepared_1_.fits'
    for p in paths:
        #first estimating the vsini from FWHM to determine whether EW will work
        spectrum= ispec.read_spectrum(p+file_name)
        template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/NARVAL.Sun.370_1048nm/template.txt.gz") #for template matching
        models, ccf = ispec.cross_correlate_with_template(spectrum, template, \
                                            lower_velocity_limit=-200, upper_velocity_limit=200, \
                                            velocity_step=1.0, fourier=False)
        c = 299792458.0 # m/s
        fwhmm = models[0].sig() * (2*np.sqrt(2*np.log(2)))  
        
        vsini = 0.662557*fwhmm - 6.119825 #from calibration - see Freckelton et al 2024
        if vsini > 10:
            print('vsini estimated to be too high for Equivalent Widths - saving estimates based on spectral type instead. Continue to synthesis :) ')
            if mask =='G2':
                t = 5700
                l = 4.0
                feh = 0.05
                a = 0.00
                vmic =ispec.estimate_vmic(t, l,feh)
            
            if mask == 'K5':
                t=4440
                l=4.6
                feh = 0.5
                a = 0.00
                vmic = ispec.estimate_vmic(t, l,feh)
            if mask == 'F3':
                t=6750
                l=4.2
                feh = 0.00
                a = 0.00
                vmic = ispec.estimate_vmic(t, l,feh)
            
            params = ['teff','logg','MH','alpha','vmic']
            values = [t,l,feh,a,vmic]
            errs = [0,0,0,0,0]
            df = pd.DataFrame(zip(params,values), columns=['Parameter','Value'])
            print(df)
            df['Errors'] = errs
        
            df.to_csv(p +'/params.csv', index=False)
        
        else:
         ####=============fitting line regions===========================

            if resolution <= 50000:
                line_regions_with_atomic_data = ispec.read_line_regions(ispec_dir + "input/regions/47000_GES/{}_ew_ispec_good_for_params_all_extended.txt".format(code))
            elif 50000< resolution <= 90000:
                line_regions_with_atomic_data = ispec.read_line_regions('solar/SOPHIE/solar_line_regions_all_ew.txt')      
            elif resolution > 90000:
                line_regions_with_atomic_data = ispec.read_line_regions('solar/HARPS/solar_line_regions_all_ew_12.txt')
        



            line_regions_with_atomic_data = line_regions_with_atomic_data[np.logical_or(line_regions_with_atomic_data['element'] == "Fe 1", line_regions_with_atomic_data['element'] == "Fe 2")]

            smoothed_star_spectrum = ispec.convolve_spectrum(spectrum, 2*to_resolution)
            line_regions_with_atomic_data = ispec.adjust_linemasks(smoothed_star_spectrum, line_regions_with_atomic_data, max_margin=0.5)
            star_continuum_model = ispec.fit_continuum(spectrum, fixed_value=1.0, model="Fixed value")
            #--- Fit the lines but do NOT cross-match with any atomic linelist since they already have that information
            linemasks = ispec.fit_lines(line_regions_with_atomic_data, spectrum , star_continuum_model,\
                                            atomic_linelist = None, \
                                            max_atomic_wave_diff = 0.005,\
                                            check_derivatives = False, \
                                            discard_gaussian=False, \
                                            smoothed_spectrum=None, \
                                            discard_voigt=True, \
                                            free_mu=True, crossmatch_with_mu=False, closest_match=False)

            # Exclude lines that have not been successfully cross matched with the atomic data
            # because we cannot calculate the chemical abundance (it will crash the corresponding routines)
            rejected_by_atomic_line_not_found = (linemasks['wave_nm'] == 0)
            linemasks = linemasks[~rejected_by_atomic_line_not_found]

            ###=============discarding any bad masks===========================

            flux_peak = spectrum['flux'][linemasks['peak']]
            flux_base = spectrum['flux'][linemasks['base']]
            flux_top = spectrum['flux'][linemasks['top']]
            bad_mask = np.logical_or(linemasks['wave_peak'] <= linemasks['wave_base'], linemasks['wave_peak'] >= linemasks['wave_top'])
            bad_mask = np.logical_or(bad_mask, flux_peak >= flux_base)
            bad_mask = np.logical_or(bad_mask, flux_peak >= flux_top)
            linemasks = linemasks[~bad_mask]

            #================Exclude lines with EW equal to zero=========

            rejected_by_zero_ew = (linemasks['ew'] == 0)
            linemasks = linemasks[~rejected_by_zero_ew]

            #================Exclude lines that may be affected by tellurics========

            rejected_by_telluric_line = (linemasks['telluric_wave_peak'] != 0)
            linemasks = linemasks[~rejected_by_telluric_line]

            #================Model spectra from EW===========================
            if mask =='G2':
                initial_teff = 5700
                initial_logg = 4.0
                initial_MH = 0.05
                initial_alpha = 0.00
                initial_vmic =ispec.estimate_vmic(initial_teff, initial_logg, initial_MH)
            if mask == 'K5':
                initial_teff=4440
                initial_logg=4.6
                initial_MH = 0.5
                initial_alpha = 0.00
                initial_vmic = ispec.estimate_vmic(initial_teff, initial_logg, initial_MH)
            if mask == 'F3':
                initial_teff=6750
                initial_logg=4.2
                initial_MH = 0.00
                initial_alpha = 0.00
                initial_vmic = ispec.estimate_vmic(initial_teff, initial_logg, initial_MH)
            max_iterations = 15

        
            model = ispec_dir + "input/atmospheres/ATLAS9.Castelli/"

            atomic_linelist_file= ispec_dir +"/input/linelists/transitions/SPECTRUM.300_1100nm/atomic_lines.tsv"

            solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
            # Load model atmospheres
            modeled_layers_pack = ispec.load_modeled_layers_pack(model)

            # Load SPECTRUM abundances
            solar_abundances = ispec.read_solar_abundances(solar_abundances_file)


            # Validate parameters
            if not ispec.valid_atmosphere_target(modeled_layers_pack, {'teff':initial_teff, 'logg':initial_logg, 'MH':initial_MH, 'alpha':initial_alpha}):
                msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] \
                    fall out of theatmospheric models."
                print(msg)

            # Reduced equivalent width
            # Filter too weak/strong lines
            # * Criteria presented in paper of GALA

            efilter = np.logical_and(linemasks['ewr'] >= -6.0, linemasks['ewr'] <= -4.3)
            # Filter high excitation potential lines
            # * Criteria from Eric J. Bubar "Equivalent Width Abundance Analysis In Moog"
            efilter = np.logical_and(efilter, linemasks['lower_state_eV'] <= 5.0)
            efilter = np.logical_and(efilter, linemasks['lower_state_eV'] >= 0.5)
            ## Filter also bad fits
            efilter = np.logical_and(efilter, linemasks['rms'] < 1.00)
            # no flux
            noflux = spectrum['flux'][linemasks['peak']] < 1.0e-10
            efilter = np.logical_and(efilter, np.logical_not(noflux))
            unfitted = linemasks['fwhm'] == 0
            efilter = np.logical_and(efilter, np.logical_not(unfitted))

            results = ispec.model_spectrum_from_ew(linemasks[efilter], modeled_layers_pack, \
                            solar_abundances, initial_teff, initial_logg, initial_MH, initial_alpha, initial_vmic, \
                            free_params=["teff", "logg","vmic"], \
                            adjust_model_metalicity=True, \
                            max_iterations=max_iterations, \
                            enhance_abundances=True, \
                            #outliers_detection = "robust", \
                            #outliers_weight_limit = 0.90, \
                            outliers_detection = "sigma_clipping", \
                            #sigma_level = 3, \
                            tmp_dir = None, \
                            code=code)

            params, errors, status, x_over_h, selected_x_over_h, fitted_lines_params, used_linemasks = results
            data = []
        
            for key,value in params.items():
                data.append([key,value])
            errs=[]
            for key,value in errors.items():
                errs.append(value)



            df = pd.DataFrame(data, columns=['Parameter','Value'])
            df['Errors']=errs
            df.to_csv(p +'/params.csv', index=False)

            stat=[]
            for  key,value in status.items():
                stat.append([key,value])
            df=pd.DataFrame(stat)
            df.to_csv(p +'/status.csv', index=False)

            df =pd.DataFrame(fitted_lines_params)
            df.to_csv(p +'/fitted_line_params.csv', index=False)

            df =pd.DataFrame(used_linemasks)
            df.to_csv(p +'/used_linemasks.csv', index=False)

            np.savetxt(p + '/x_over_h_ew.csv',x_over_h, delimiter=',')
            file=open(p +'/selected_x_over_h.txt', 'w')
            for element in selected_x_over_h:
                file.write(str(element)+ '\n')
            x_h = [i for i in x_over_h if str(i) != 'nan']
            lower_state_e=[row[4] for row in used_linemasks]
            element = [row[0] for row in used_linemasks]
            lower_state_ev_1 = [lower_state_e[i] for i in range(len(lower_state_e)) if element[i]=='Fe 1'] 
            lower_state_ev_2 = [lower_state_e[i] for i in range(len(lower_state_e)) if element[i]=='Fe 2' ] 

        
        
            avg = [mean(x_h) for i in range(len(lower_state_e))]
            x_h_1 = [x_h[i] for i in range(len(x_h)) if element[i]=='Fe 1' ] 

            x_h_2 = [x_h[i] for i in range(len(x_h)) if element[i]=='Fe 2' ] 


    return 
run = widgets.Button(description='Run Analysis Preparation')
out9 = widgets.Output(layout={'border': '1px solid black'})
box9 = widgets.VBox([widgets.VBox([run]),out9])

run.on_click(parameters_from_ew)

display(box9)
