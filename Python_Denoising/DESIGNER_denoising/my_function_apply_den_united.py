#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 24/04/2023

@author: Maxime Yon

"""

import nibabel as ni 
import numpy as np
import os
import sys
import shutil
import denoise as mpden

# my_home_path = "C:/Users/Administrateur/Mon_Drive/Matlab/ProcessingPV360/data/"
my_data_path = str(sys.argv[1]);


if not os.path.exists(my_data_path + "/orig"):
    os.makedirs(my_data_path + "/orig")
    shutil.move(my_data_path + "/data.nii.gz", my_data_path + "/orig/data.nii.gz")

nii = ni.load(my_data_path + "/orig/data.nii.gz")
dwi = np.array(nii.dataobj)

# for 2D
#Signal, Sigma,Nparameters = mpden.denoise(dwi, kernel=[15,15,1],patchsize=125,shrinkage='frobenius')
# for 3D
Signal, Sigma,Nparameters = mpden.denoise(dwi, kernel=[7,7,7],patchsize=125,shrinkage='frobenius')

niout = ni.Nifti1Image(Signal, nii.affine, nii.header)
niout.header.set_data_dtype(np.float32)

if not os.path.exists(my_data_path + "/denoised"):
    os.makedirs(my_data_path + "/denoised")

ni.save(niout, my_data_path + "/denoised/data.nii.gz")
nii.header["dim"][0] = 3
nii.header["dim"][4] = 1