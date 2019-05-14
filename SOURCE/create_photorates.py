#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
file name: create_photorates.
author: Sacha Gavino
date: May 2019
language: PYTHON 3.7
version of program: 1.0
main: create a table of photorate values for all species generated by the ISRF. The photorate values are given before being weighted VERTICALLY by the disk medium.
comments: This code is meant to be implemented in disk_structure.py
"""
import os
import sys
import glob

import numpy as np
import pandas as pd
import re


#-- FUNCTION 1: create a list of species computed.
def X():
    files = [f for f in os.listdir('cross-sections_diss/') if re.findall('.txt', f)]
    X = [re.sub('.txt', '', fname) for fname in files]

    return X
#------------------------------------------------

#_________________________________________________#
#                      MAIN                       #
#_________________________________________________#
if __name__ == "__main__":

    #---------variables---------
    field = 'standard_DRAINE_ISRF_1978' 
    X = X()
    list_photorates = []
    head="species\tphotorate(s-1)"
    #---------------------------

    for species in X: 

        #-----read cross-sections tables-----
        cross_sections = pd.read_table('cross-sections_diss/%s.txt' %  species, sep=" ", comment='#', header=None)
        name_cs = ['wv', 'abs', 'diss', 'ion'] # column names: wavelength(wv), absorption(abs), dissociation(diss), ionization(ion)
        cross_sections.columns = name_cs
        nb_rows=len(cross_sections)
        #------------------------------------

        #-----read UV flux tables-----
        flux_isrf = pd.read_table('radiation_fields/%s.txt' %  field, sep=" ", comment='#', header=None)
        name_flux = ['wv', 'flux'] # column names: wavelength(wv), UV flux(flux)
        flux_isrf.columns = name_flux
        #-----------------------------

        #-----create the step of integration dlambda-----
        dlambda = (flux_isrf['wv'][11] - flux_isrf['wv'][10]) #Since dlambda is constant we use arbitrary rows to compute it (here 11 and 10).
        #------------------------------------------------

        #-----Make sure we have the same size and first value of wavelength in both flux and cross-sections-----
        start = cross_sections['wv'].iloc[0]
        #print(start)
        flux_isrf = flux_isrf[flux_isrf.wv >= start]
        #-------------------------------------------------------------------------------------------------------

        #-----multiply cross-section table with the flux and sum-----
        product=cross_sections['diss']*flux_isrf['flux']*dlambda
        product_filtered = product.dropna()                     ####!!!! we will get the wavelength at which the cs table starts, then we cut everything below that wv in the flux table. easy
        summed = product.sum()
        list_photorates.append(summed)
        #------------------------------------------------------------

    #-----create photorates.in file-----
    photorates = pd.Series(list_photorates)
    spe = pd.Series(X)
    photorates_table = pd.concat([spe, photorates], axis=1)
    #print(photorates_table.values)
    np.savetxt('photorates.in', photorates_table.values, fmt='%-11s %.6e', newline='\n', header=head, comments='#', encoding=None)    
    #-----------------------------------
