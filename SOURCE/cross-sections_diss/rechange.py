#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
file name: create_photorates.
author: Sacha Gavino
date: May 2019
language: PYTHON 3.7
version of program: 1.0
main: rechange resolution to fit with the uv flux data.
"""
import os
import sys
import glob

import numpy as np
import pandas as pd
import re

name = ['wv', 'abs', 'diss', 'ion']

for filename in os.listdir():
    if filename.endswith("notres.txt"):
        adapted=pd.DataFrame()
        X = re.sub('_notres.txt', '', filename)
        cross_sections = pd.read_table(filename, sep=" ", comment='#', header=None)
        cross_sections.columns = name

        adapted = cross_sections.iloc[::10, :]
        
        header='%s cross sections\n\
#photoabsorption   -- Photoabsorption cross section (cm2)\n\
#photodissociation -- Photodissociation cross section (cm2)\n\
#photoionisation   -- Photoionisation cross section (cm2)\n\
#This version has been resampled to a uniform 0.1nm grid while maintaining a constant integrated cross section (using trapezoidal integration).\n\
#wavelength photoabsorption photodissociation photoionisation' % X

        np.savetxt('%s.txt' % X, adapted, fmt='%.6e', delimiter=' ', newline='\n', header=header , comments='# ', encoding=None)
        