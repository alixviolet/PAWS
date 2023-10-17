# spectral_analysis
	If using linux/macOS, install iSpec directly from the guides at https://www.blancocuaresma.com/s/iSpec/manual/installation/linux or https://www.blancocuaresma.com/s/iSpec/manual/installation/osx/anacondapython.
 	Make sure all of the packages from the list at https://www.blancocuaresma.com/s/iSpec/manual/installation/linux/apt are installed.
 
	If using windows, then you will need to set up a virtual machine to perform analysis. I did this using https://www.virtualbox.org/ and following the guide on https://linuxconfig.org/how-to-install-ubuntu-20-04-on-virtualbox. If the window refuses to go into fullscreen mode, follow the guide on https://www.itzgeek.com/post/how-to-install-virtualbox-guest-additions-on-ubuntu-20-04/.
	To install iSpec in the virtual machine, follow the guide at https://www.blancocuaresma.com/s/iSpec/manual/installation/linux/anacondapython. Make sure all of the packages from the list at https://www.blancocuaresma.com/s/iSpec/manual/installation/linux/apt are installed. 
	All functions necessary for spectral analysis are contained within the .py files, with the PAWS jupyter notebook set up as a GUI to utilise these.
	It is important to check the format of your spectra - the code has been optimised for use with SOPHIE spectra. Depending on the instrument, you may be working with spectra of a different resolution so make sure to check this (for both HE and HR spectra if relevant). Also, some spectra may already be coadded, normalised or rv/telluric/barycentric velocity corrected - if so, you can select yes/no as appropriate in the GUI.

**Using the jupyter GUI**

First, run the 'import setup cell'.
Add the paths to iSpec and your data, and hit the corresponding buttons.
Select the appropriate options for whether you are using single targets, their spectral type, and whether coadding is required.

Input the resolution of the observations, and your minimum & maximum wavelengths to analyse between.
Make sure your individual spectra have a common string, so the code can single these out from other files.
Select the correct wavelength unit for your spectra.

Run the 'import coadding' cell & press the button to begin the normalisation & rv correction process, which will then coadd the spectra if necessary.

run 'import pipeline' & press the button to perform EW analysis, which will save a list of parameters for the star.

Run 'import synth' & press the button to perform spectral synthesis, which requires EW analysis to be run prior, as it uses the parameters as inputs. 
