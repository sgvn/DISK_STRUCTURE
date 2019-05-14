#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
file name: create_photorates.
author: Sacha Gavino
date: May 2019
language: PYTHON 3.7
version of program: 1.0
main: reround the files which rounding didn't work. can be deleted afterwards.
"""
import os
import sys
import glob

import numpy as np
import pandas as pd
import re

X = 'HCO+'
name = ['wv', 'abs', 'diss', 'ion'] # column names

# import cross_sections data
cross_sections = pd.read_table('%s_0.1nm.txt' % X, sep=" ", comment='#', header=None) # 2)
cross_sections.columns = name

wv = round(cross_sections['wv'], 1)


abso = cross_sections['abs']
diss = cross_sections['diss']
ion = cross_sections['ion']
rounded = np.stack((wv, abso, diss, ion), axis=-1)
print(rounded)
header='%s cross sections\n\
photoabsorption   -- Photoabsorption cross section (cm2)\n\
photodissociation -- Photodissociation cross section (cm2)\n\
photoionisation   -- Photoionisation cross section (cm2)\n\
This version has been resampled to a uniform 0.1nm grid while maintaining a constant integrated cross section (using trapezoidal integration).\n\
wavelength photoabsorption photodissociation photoionisation' % X
np.savetxt('%s_notres.txt' % X, rounded, fmt='%.6e', delimiter=' ', newline='\n', header=header , comments='# ', encoding=None)