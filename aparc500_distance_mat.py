#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to generate distance matrices for BrainSMASH with 500aparc parcellation
@author: C. Schleifer
"""

#import pandas as pd
#import numpy as np
#pd.set_option('display.max_rows', 1000)
#pd.set_option('display.max_columns', 500)

# install brainsmash
# pip install brainsmash
# import relevant functions
#from brainsmash.mapgen.base import Base
#from brainsmash.utils.dataio import load
from brainsmash.workbench.geo import cortex
#from brainsmash.workbench.geo import subcortex
#from brainsmash.workbench.geo import parcellate

# compute cortical distance matrices for each hemisphere
# https://brainsmash.readthedocs.io/en/latest/gettingstarted.html#computing-a-cortical-distance-matrix    print(hemi)
label="lh.500.aparc.dlabel.nii"
surface = "S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii"
cortex(surface=surface, dlabel=label, outfile="lh.500.aparc.distmat.txt", euclid=False)
