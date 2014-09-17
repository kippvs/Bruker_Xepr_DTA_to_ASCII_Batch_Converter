# Bruker Xepr DTA to ASCII Batch Converter  
Filename:       BrukerASCIIConvert.py  
Author:         Kipp J van Schooten (kippvs@gmail.com)  
Version:        0.1 beta  
Date:           2014-07-19 - 2014-07-30  

## Description:  
I've worked with a Bruker Elexsys E580 Electron Spin Resonance machine for a 
few years now on a daily basis and have acquired thousands of data sets. Besides 
being quite 'crashy,' the control software, Xepr, does not allow the user to 
perform a batch conversion of data sets to ASCII format for external analysis. 
Doing this one file at a time regularly added 15-20 min. to 
each measurement day. This Python script, and associated IPython Notebook, are 
meant to avoid this particular frustration and keep this chore down to about 20 
seconds.

## Usage:  
* user pastes file directory to recurse through
* matches to .DTA and .DSC are found
* list of matches is itereated through
	- binary (BE decimal) data is loaded from DTA file as double
	- sweep params are loaded from DSC file
	- composite dataset created
	- export to .txt file performed

## Notes:  
* data files must have BOTH a .DTA and .DSC file in same directory
* will overwrite existing .txt files of same name
* supports 1D, 2D, and 3D data sets, or either real or complex format
* all output data columns populated with floats
* execution time limited by HDD access speeds (i.e: fileIO)
* corrupt DTA files usually result in ValueError (mismatched axis arrays)
* offered AS IS, but feel free to contact me with questions/bugs
