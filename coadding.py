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



from setup import datapathWidget,number,coadd, res_widget,ispecpathWidget, max_widget, min_widget, s1d_widget,unit
ispec_dir = ispecpathWidget.value

if ispec_dir[-1] != '/':
            ispec_dir=ispec_dir+'/'
import ispec 
c_str=s1d_widget.value

def degrade_norm(b):
    print('normalising')
    resolution=res_widget.value
    #path=[]
    wavelength_min = min_widget.value
    wavelength_max=max_widget.value
    
    if number.value=='Multiple' and coadd.value=='No':
        paths = [ f.path for f in os.scandir(datapathWidget.value) if f.is_dir()]
    if number.value=='Single' and coadd.value=='No':
        paths = [datapathWidget.value]
    print(paths)
    with out8:
        for p in paths:
            #print([f.name for f in os.scandir(p+'/') if f.is_file()])
            raw_spectra = [f.name for f in os.scandir(p+'/') if f.is_file() and str(c_str) in f.name] #and 'S1D' in f.name]
            print(raw_spectra)
            cor_spectrum = ispec.read_spectrum(p+'/'+raw_spectra[0])
            h=fits.open(p+'/'+raw_spectra[0])
            print(h)

            print(cor_spectrum['flux'])
            #cor_spectrum['flux'] = [i/1e-13 for i in cor_spectrum['flux']]
            print(cor_spectrum['flux'])
            plt.plot(cor_spectrum['waveobs'],cor_spectrum['flux'])
            if str(unit.value) == 'Angstroms':
                cor_spectrum['waveobs']=cor_spectrum['waveobs']/10
            snr = ispec.estimate_snr(cor_spectrum['flux'], num_points=20)  
            cor_spectrum['err'] = cor_spectrum['flux']/snr
                                
            logging.info("Radial velocity determination with template...")
            template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/NARVAL.Sun.370_1048nm/template.txt.gz")
            models, ccf = ispec.cross_correlate_with_template(cor_spectrum, template, \
                                            lower_velocity_limit=-1000, upper_velocity_limit=1000, \
                                            velocity_step=1.0, fourier=False)
            components = len(models)
            print(components)
                    # First component:
            rv = np.round(models[0].mu(), 2) # km/s
            rv_err = np.round(models[0].emu(), 2) # km/s
            cor_spectrum = ispec.correct_velocity(cor_spectrum, rv)
            #if resolution > 47000:
             #   to_resolution = 47000
              #  deg_spectrum = ispec.convolve_spectrum(cor_spectrum, to_resolution, \
                                                                #from_resolution=resolution)
               # print('degraded resolution')
            #else:
                #to_resultion = resolution
                #deg_spectrum=cor_spectrum
                #print('kept resolution the same')
            deg_spectrum = cor_spectrum
          #  deg_spectrum=cor_spectrum
                                     #--- ontinuum fit 
            model = "Splines" # "Polynomy"
            degree = 2
            nknots = None # Automatic: 1 spline every 5 nm
            from_resolution = to_resolution = res_widget.value

            # Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
            order='median+max'
            median_wave_range=0.05
            max_wave_range=1.0
            plt.scatter(deg_spectrum['waveobs'], deg_spectrum['flux'])
            plt.show()
            star_continuum_model = ispec.fit_continuum(deg_spectrum, from_resolution=from_resolution, \
                                                        nknots=nknots, degree=degree, \
                                                        median_wave_range=median_wave_range, \
                                                        max_wave_range=max_wave_range, \
                                                        model=model, order=order, \
                                                        automatic_strong_line_detection=True, \
                                                        strong_line_probability=0.5, \
                                                        use_errors_for_fitting=True)
                            #--- Normalize -------------------------------------------------------------
            normalized_star_spectrum = ispec.normalize_spectrum(deg_spectrum, star_continuum_model,     consider_continuum_errors=False)
            
            print('done normalising')
            file_name=p +'/prepared_1_.fits'
            print(file_name)
            ispec.write_spectrum(normalized_star_spectrum, file_name)
            #print('writing')
go = widgets.Button(description='Prepare Spectra')
go.style.button_color = '#ffbc3c'
out7 = widgets.Output(layout={'border': '1px solid black'})
box7 = widgets.VBox([widgets.VBox([go]),out7])


nom = widgets.Button(description='Normalise Spectra')
nom.style.button_color = '#ffbc3c'
out8 = widgets.Output(layout={'border': '1px solid black'})
box8 = widgets.VBox([widgets.VBox([nom]),out8])

nom.on_click(degrade_norm)

if coadd.value=='Yes':
    display(box7)
elif coadd.value=='No':
    display(box8)
def coadd_spectra(b):
    print('Coadding...')
    resolution=res_widget.value
    paths=[]
    to_resolution=resolution #only for now
    wavelength_min = min_widget.value
    wavelength_max=max_widget.value
    if number.value=='Multiple' and coadd.value=='Yes':
        paths = [f.path for f in os.scandir(datapathWidget.value) if f.is_dir()]
        #print(paths)
    if number.value=='Single' and coadd.value=='Yes':
        paths = [datapathWidget.value]
    with out7:
        for p in paths:
            raw_spectra = [f.name for f in os.scandir(p) if f.is_file() and str(c_str) in f.name]
            #print(raw_spectra)
            spectra=[]
            snr_list=[]
            for i in range(len(raw_spectra)):
                print(raw_spectra[i])
                spectrum = ispec.read_spectrum(p+'/'+raw_spectra[i])
                if str(unit.value) == 'Angstroms':
                    spectrum['waveobs']=spectrum['waveobs']/10
                print('read spectrum!')
                hdul = fits.open(p+'/'+raw_spectra[i])
                hdr=hdul[0].header
                #spectrum['flux'] = [i/1e-13 for i in spectrum['flux']]
                snr = ispec.estimate_snr(spectrum['flux'], num_points=10)                              
                spectrum['err'] = spectrum['flux']/snr
                snr_list.append(snr)
                if snr> 10:
                    
                    logging.info("Radial velocity determination with template...")
                    template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/NARVAL.Sun.370_1048nm/template.txt.gz")
                    models, ccf = ispec.cross_correlate_with_template(spectrum, template, \
                                            lower_velocity_limit=-200, upper_velocity_limit=200, \
                                            velocity_step=1.0, fourier=False)
                    components = len(models)
                    print(components)
                    # First component:
                    rv = np.round(models[0].mu(), 2) # km/s
                    rv_err = np.round(models[0].emu(), 2) # km/s
                    cor_spectrum = ispec.correct_velocity(spectrum, rv)

                    #if resolution > 47000:
                     #   to_resolution = 47000
                      #  deg_spectrum = ispec.convolve_spectrum(cor_spectrum, to_resolution, \
                                                               # from_resolution=resolution)
                    #else:
                     #   to_resultion = 40000
                      #  deg_spectrum=cor_spectrum
                    deg_spectrum=cor_spectrum
                    
                                 #--- Continuum fit -------------------------------------------------------------
                    model = "Splines" # "Polynomy"
                    degree = 2
                    nknots = None # Automatic: 1 spline every 5 nm
                    from_resolution = to_resolution

                    # Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
                    order='median+max'
                    median_wave_range=0.05
                    max_wave_range=1.0

                    star_continuum_model = ispec.fit_continuum(deg_spectrum, from_resolution=from_resolution, \
                                                nknots=nknots, degree=degree, \
                                                median_wave_range=median_wave_range, \
                                                max_wave_range=max_wave_range, \
                                                model=model, order=order, \
                                                automatic_strong_line_detection=True, \
                                                strong_line_probability=0.5, \
                                                use_errors_for_fitting=True)
                    #--- Normalize -------------------------------------------------------------
                    normalized_star_spectrum = ispec.normalize_spectrum(deg_spectrum, star_continuum_model,     consider_continuum_errors=False)
                    spectra.append(normalized_star_spectrum)
            wavelengths = np.arange(wavelength_min, wavelength_max, 0.001)
            flux=[]
            err=[]
            for i in spectra:
                if len(spectra)>=1:
                    resampled_star_spectrum = ispec.resample_spectrum(i, wavelengths, method='linear', zero_edges=True)

                    flux.append(resampled_star_spectrum['flux'])
                    err.append(resampled_star_spectrum['err']) 
            if len(spectra)>=1:
                coadded_spectrum = ispec.create_spectrum_structure(wavelengths, flux = np.average(flux, axis=0),err=np.average(err, axis=0))
                spec_snr = ispec.estimate_snr(coadded_spectrum['flux'], num_points=10)
                smoothed_star_spectrum = ispec.convolve_spectrum(coadded_spectrum, to_resolution)
                snr_list.append(spec_snr)
                file=open(p +'/SNR_list.txt', 'w')
                for element in snr_list:
                    file.write(str(element)+ '\n')
                file_name=p +'/prepared_1_.fits'
                ispec.write_spectrum(coadded_spectrum, file_name)
            else:
                file=open(p+'/Spectra_too_low_SNR', 'w')
                for element in snr_list:
                    file.write(str(element)+ '\n')
        return 
go.on_click(coadd_spectra)


    

    