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


from setup import datapathWidget,number,coadd, res_widget,ispecpathWidget, input_widget, max_widget, min_widget
ispec_dir=ispecpathWidget.value

if ispec_dir[-1] != '/':
            ispec_dir=ispec_dir+'/'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec 







def merged_synthesis(b):
    code='spectrum'
    
    resolution=res_widget.value
    mask = input_widget.value
    if resolution >= 47000:
        to_resolution = 47000
                
    elif resolution <47000:
        to_resolution = resolution
    if number.value=='Multiple':
        paths = [ f.path for f in os.scandir(datapathWidget.value) if f.is_dir()]
    if number.value=='Single':
        paths = [datapathWidget.value]   
    file_name = '/prepared_1_.fits'
    print(paths)
    for p in paths:
        spectrum= ispec.read_spectrum(p+file_name)


        #spectrum = ispec.read_spectrum(merged_spectrum)
        #print('read spectrum')


        #print('assigned res')
        #===========Model spectra===========================
        # Parameters
        ew_res = df=pd.read_csv(p+'/params.csv')
        initial_teff =ew_res.loc[0]['Value']
        #initial_teff=6302
        initial_logg = ew_res.loc[1]['Value']
        #initial_logg=4.260822131426127
        initial_MH = ew_res.loc[2]['Value']
        #initial_MH = 0.767009374659
        if initial_MH >0.5:
            initial_MH = 0.5
        initial_alpha = ew_res.loc[3]['Value']#ispec.determine_abundance_enchancements(initial_MH) #use input from ew
        initial_vmic = ew_res.loc[4]['Value']
        #ispec.estimate_vmic(initial_teff, initial_logg, initial_MH)
        initial_vmac = ispec.estimate_vmac(initial_teff, initial_logg, initial_MH)
        if initial_vmac>50:
            initial_vmac=50
        #ispec.estimate_vmac(initial_teff, initial_logg, initial_MH)
        initial_vsini = 2
        #ispec.estimate_vsini(initial_teff, initial_logg, initial_MH)
        #1.27332434061464
        initial_limb_darkening_coeff = 0.6
        initial_R = to_resolution
        initial_vrad = 0
        max_iterations = 6 #CHANGE THIS
        print('set params')
        # Selected model amtosphere, linelist and solar abundances

        #model = ispec_dir + "/input/atmospheres/MARCS.GES/"
        model = ispec_dir + "input/atmospheres/ATLAS9.Castelli/"
        print('model selected')

        #GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv"

        print('linelist selected')
        if "ATLAS" in model:
            solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
        else:
            # MARCS
            solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"
        print('solar abundances loaded')
        isotope_file = ispec_dir + "/input/isotopes/SPECTRUM.lst"
        atomic_linelist_file=ispec_dir +"/input/linelists/transitions/SPECTRUM.300_1100nm/atomic_lines.tsv"
        
        # Load chemical information and linelist
        atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, wave_base=np.min(spectrum['waveobs']), wave_top=np.max(spectrum['waveobs']))
        atomic_linelist = atomic_linelist[atomic_linelist['theoretical_depth'] >= 0.01] # Select lines that have some minimal contribution in the sun

        isotopes = ispec.read_isotope_data(isotope_file)
        print('linelist and isotope files read')
        # Load model atmospheres
        modeled_layers_pack = ispec.load_modeled_layers_pack(model)

        # Load SPECTRUM abundances
        solar_abundances = ispec.read_solar_abundances(solar_abundances_file)
        print('model and abundances loaded')

        # Free parameters - these can be changed
        #free_params = ["teff", "logg", "MH", "vmic", "vmac", "vsini", "R", "vrad", "limb_darkening_coeff"]



        free_params = ["vsini","teff", "logg"]
        print('free parameters set')
        # Free individual element abundance

        chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
        chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
        #free_abundances = ispec.create_free_abundances_structure(["He"], chemical_elements, solar_abundances)
        #free_abundances['Abund'] += initial_MH # Scale to metallicity
        free_abundances = None
        linelist_free_loggf = None

        # Line regions
        line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_GES/{}_synth_good_for_params_all_extended.txt".format(code))
        print('line regions read')
        ## Select only some lines to speed up the execution if desired, but this isn't recommended
       # line_regions = line_regions[np.logical_or(line_regions['note'] == 'Ti 1', line_regions['note'] == 'Ti 2')]
        #line_regions = ispec.adjust_linemasks(spectrum, line_regions, max_margin=0.5)
        #line_regions = line_regions[np.logical_or(line_regions['note'] == 'Fe 1', line_regions['note'] == 'Fe 2')]
        line_regions = ispec.adjust_linemasks(spectrum, line_regions, max_margin=0.5)
        # Read segments if we have them or...
        #segments = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
        # ... or create the segments
        segments = ispec.create_segments_around_lines(line_regions, margin=0.25)

        print('segments created')

        star_continuum_model = ispec.fit_continuum(spectrum, fixed_value=1.0, model="Fixed value")
        print('continuum fit')

        obs_spec, modeled_synth_spectrum, params, errors, abundances_found, loggf_found, status, stats_linemasks = \
                ispec.model_spectrum(spectrum, star_continuum_model, \
                modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, free_abundances, linelist_free_loggf, initial_teff, \
                initial_logg, initial_MH, initial_alpha, initial_vmic, initial_vmac, initial_vsini, \
                initial_limb_darkening_coeff, initial_R, initial_vrad, free_params, segments=segments, \
                linemasks=line_regions, \
                enhance_abundances=True, \
                use_errors = True, \
                vmic_from_empirical_relation = False, \
                vmac_from_empirical_relation = True, \
                max_iterations=max_iterations, \
                tmp_dir = None, \

                code=code)

        #chisquared of how each line fits 
        print(errors)
        data = []
        for key,value in params.items():
            data.append([key,value])
        errs=[]
        for key,value in errors.items():
            errs.append(value*3)
        errs[0] = (np.sqrt(errs[0]**2 +100**2))
        print(errs)


        df = pd.DataFrame(data, columns=['Parameter','Value'])
        df['Errors']=errs
        mh = df.at[2, 'Value']
        alpha = df.at[3, 'Value']

        
        df.to_csv(p +'/params_synth_pipeline.csv', index=False)

        stat=[]
        for  key,value in status.items():
            stat.append([key,value])
        df=pd.DataFrame(stat)
        df.to_csv(p +'/status_synth_pipeline.csv', index=False)

        df =pd.DataFrame(stats_linemasks)
        df.to_csv(p +'/stats_linemasks_synth_pipeline.csv', index=False)


        wavelengths = np.arange(min_widget.value, max_widget.value, 0.001)
        resampled_star_spectrum = ispec.resample_spectrum(modeled_synth_spectrum, wavelengths, method='linear', zero_edges=True)
        file_name_save=p +'/synthesised_spectrum_pipeline.fits'
        ispec.write_spectrum(resampled_star_spectrum, file_name_save)

        print(params, errors)
    return
runn = widgets.Button(description='Run Analysis')
outn = widgets.Output(layout={'border': '1px solid black'})
boxn = widgets.VBox([widgets.VBox([runn]),outn])

runn.on_click(merged_synthesis)

display(boxn)
