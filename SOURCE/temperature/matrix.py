#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os

radius="040"
size=15
spatial_pt=63
temperature = np.loadtxt('T%s.txt' % radius, delimiter=' ', unpack=True)
print(temperature[size,spatial_pt])  #column, row