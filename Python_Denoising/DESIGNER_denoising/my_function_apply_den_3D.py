#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 24/04/2023

@author: Maxime Yon

"""

import nibabel as ni 
import numpy as np
import sys
import denoise as mpden

# my_home_path = "C:/Users/Administrateur/Mon_Drive/Matlab/ProcessingPV360/data/"
my_data_path = str(sys.argv[1]);
# my_data_path = "ON-81-mbti-pilot-3d/17"

nii = ni.load(my_data_path + "/dataRaw.nii.gz")
dwi = np.array(nii.dataobj)
nrep = dwi.shape[3]; "dimensions start from 0"

Signal, Sigma,Nparameters = mpden.denoise(dwi, kernel=[7,7,7],patchsize=125,shrinkage='frobenius')

niout = ni.Nifti1Image(Signal, nii.affine, nii.header)
niout.header.set_data_dtype(np.float32)

ni.save(niout,  my_data_path + "/dataDen.nii.gz")
nii.header["dim"][0] = 3
nii.header["dim"][4] = 1

Sigmaout = ni.Nifti1Image(Sigma, nii.affine, nii.header)
Sigmaout.header.set_data_dtype(np.float32)

ni.save(Sigmaout,  my_data_path + "/sigma.nii.gz")
nii.header["dim"][0] = 3
nii.header["dim"][4] = 1
