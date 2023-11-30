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


#-----Getting iSpec Path and Importing-----------
ispecpathWidget = widgets.Text(
    placeholder='Enter iSpec Path',
    description='iSpec Path:',
    disabled=False,
    style={'description_width': 'initial'},
    button_style=''
)

get_path=widgets.Button(description='Import iSpec')
get_path.style.button_color = '#ffbc3c'
out3=widgets.Output(layout={'border': '1px solid black'})
box3=widgets.HBox([widgets.HBox([ispecpathWidget, get_path]), out3])


def set_ispec(b):
    with out3:
        ispec_dir = ispecpathWidget.value
        if ispec_dir[-1] != '/':
            ispec_dir=ispec_dir+'/'
        sys.path.insert(0, ispec_dir)
        import ispec 
        print('iSpec Imported')
        return

get_path.on_click(set_ispec)
#data_path = get_path.on_click(set_path)
display(box3)
#-----Getting Data Path--------------
datapathWidget=widgets.Text(
    placeholder='Enter Data Path',
    description='Data Path:',
    disabled=False,
    style={'description_width': 'initial'}
)
get_path = widgets.Button(description='Select Data Path')
get_path.style.button_color = '#ffbc3c'
outt = widgets.Output(layout={'border': '1px solid black'})
boxx=widgets.HBox([widgets.HBox([datapathWidget, get_path]), outt])
def set_path(c):
    global path
    path = datapathWidget.value
    with outt:
        path = datapathWidget.value
        print('Data Path Set')
    return path
get_path.on_click(set_path)
display(boxx)


p= datapathWidget.value
print(p)


#---------Multiple or one spectrum?-------

print('Are you analysing one target or multiple targets?')
number=widgets.ToggleButtons(
    options=['Multiple', 'Single'],
    description='',
    disabled=False,
    button_style='danger', 
    default=None,# 'success', 'info', 'warning', 'danger' or ''
    #tooltips=['Description of slow', 'Description of regular', 'Description of fast'],
     #icons=['check'] * 3
)
select=widgets.Button(description='Done')
select.style.button_color = '#ffbc3c'
out=widgets.Output(layout={'border': '1px solid black'})

def selection(a):
    global num
    num = str(number.value)
    with out:
        print('You have selected: '+str(number.value))
        if str(number.value) == 'Multiple':
            print('Ensure all spectra are within individual folders in your data path.')
        elif str(number.value)== 'Single':
            print('Ensure all data is within your data path.')
        
    return 
            
    
    
select.on_click(selection)
box=widgets.VBox([widgets.VBox([number, select]),out])
display(box)
#print(num)

#----------------input parameters----------------
print('Which spectral type would you like to base input parameters upon? This will be used for all targets.')
input_widget = widgets.ToggleButtons(
    options=['G2','K5', 'F3'],
    description='',
    disabled=False,
    button_style='danger',
    default=None,
)

mask = widgets.Button(description='Done')
mask.style.button_color = '#ffbc3c'
out_mask = widgets.Output(layout={'border': '1px solid black'})
box_mask=widgets.VBox([widgets.VBox([input_widget, mask]),out_mask])
def masks(b):
    with out_mask:
        print('You have selected the ' +input_widget.value +' mask.')
    return

mask.on_click(masks)
display(box_mask)
#--------------Is Coadding Necessary?----------
print('Do your spectra require coadding?')
coadd=widgets.ToggleButtons(
    options=['Yes', 'No'],
    description='',
    disabled=False,
    button_style='danger', 
    default=None,# 'success', 'info', 'warning', 'danger' or ''
    #tooltips=['Description of slow', 'Description of regular', 'Description of fast'],
     #icons=['check'] * 3
)
coadd.style.button_color = '#ffbc3c'
coaddd=widgets.Button(description='Done')
coaddd.style.button_color = '#ffbc3c'

out1=widgets.Output(layout={'border': '1px solid black'})
box1=widgets.VBox([widgets.VBox([coadd, coaddd]),out1])
def selectadd(a):
    with out1:
        print('You have selected: '+str(coadd.value))
        if str(coadd.value)== 'No':
            print('Coadding will not be performed on your spectra.')
            
    
coaddd.on_click(selectadd)
display(box1)

print('What is the spectral resolution?')
res_widget=widgets.IntText(
    value=None,
    description='Resolution:',
    disabled=False
)
sel_res = widgets.Button(description='Set Resolution')
sel_res.style.button_color = '#ffbc3c'
out4 = widgets.Output(layout={'border': '1px solid black'})
box4 = widgets.HBox([widgets.HBox([res_widget, sel_res]),out4])
def set_res(a):
    with out4:
        print('Your resolution is set to :' +str(res_widget.value))
    return
sel_res.on_click(set_res)
display(box4)
print('What are your minimum and maximum wavelengths (nm)?')
min_widget=widgets.IntText(
    
    description='Minimum :',
    disabled=False
)
max_widget=widgets.IntText(
    
    description='Maximum :',
    disabled=False
)
set_wav= widgets.Button(description='Done')
set_wav.style.button_color = '#ffbc3c'
out6 = widgets.Output(layout={'border': '1px solid black'})
box6 = widgets.HBox([widgets.HBox([min_widget, max_widget, set_wav]),out6])
def waves(a):
    print('Wavelength Range Selected.')
set_wav.on_click(waves)
display(box6)

print('What is the common string in your spectrum files? e.g. \'S1D\'')
s1d_widget=widgets.Text(
    value=None,
    description='String:',
    disabled=False
)
sel_str = widgets.Button(description='Set String')
sel_str.style.button_color = '#ffbc3c'
out7 = widgets.Output(layout={'border': '1px solid black'})
box7 = widgets.HBox([widgets.HBox([s1d_widget, sel_str]),out7])
def set_str(a):
    with out7:
        print('Your common string is :' +str(s1d_widget.value))
    return
sel_str.on_click(set_str)
display(box7)


print('What is the wavelength unit of your spectra?')
unit=widgets.ToggleButtons(
    options=['Nanometres', 'Angstroms'],
    description='',
    disabled=False,
    button_style='danger', 
    default=None,
)
#'coadd.style.button_color = '#ffbc3c'
unitt=widgets.Button(description='Done')
unitt.style.button_color = '#ffbc3c'

out8=widgets.Output(layout={'border': '1px solid black'})
box8=widgets.VBox([widgets.VBox([unit, unitt]),out8])
def sel_unit(a):
    with out8:
        print('You have selected: '+str(unit.value))

            
    
unitt.on_click(sel_unit)
display(box8)


#OS to determine width or moog

print('What is your OS?')
osx=widgets.ToggleButtons(
    options=['Ubuntu/Linux', 'Mac/iOS'],
    description='',
    disabled=False,
    button_style='danger', 
    default=None,
)

osxx=widgets.Button(description='Done')
osxx.style.button_color = '#ffbc3c'

out9=widgets.Output(layout={'border': '1px solid black'})
box9=widgets.VBox([widgets.VBox([osx, osxx]),out9])
def sel_os(a):
    with out9:
        if osx.value == 'Ubuntu/Linux':
            print('You have selected: '+str(osx.value) +', the width radiative transfer code will be used.')
        else:
            print('You have selected: ' + str(osx.value)+', the moog radiative transfer code will be used.')

            
    
osxx.on_click(sel_os)
display(box9)


    
    