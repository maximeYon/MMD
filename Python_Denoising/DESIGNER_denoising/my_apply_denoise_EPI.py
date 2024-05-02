#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 24/04/2023

@author: Maxime Yon

"""

import nibabel as ni 
import numpy as np
import os
import denoise as mpden

# my_home_path = "C:/Users/Administrateur/Mon_Drive/Matlab/ProcessingPV360/data/"
my_home_path = "C:/Users/User/Mon_Drive/Matlab/ProcessingPV360/data/"
my_data_path =  "EPI_RARE_MSME/16"

nii = ni.load(my_home_path + my_data_path + "/pdata_mdd/nii_xps/data.nii.gz")
dwi = np.array(nii.dataobj)
nrep = dwi.shape[3]; "dimensions start from 0"

DwiDen, Sigma,Nparameters = mpden.denoise(dwi, kernel=[15,15,5],patchsize=125,shrinkage='frobenius')

#DwiDen, Sigma,Nparameters = mpden.denoise(dwi, kernel=[7,7,7],patchsize=125,shrinkage='frobenius')

niout = ni.Nifti1Image(DwiDen, nii.affine, nii.header)
niout.header.set_data_dtype(np.float32)


if not os.path.exists(my_home_path + my_data_path + "/pdata_mdd/nii_xps/denoised"):
    os.makedirs(my_home_path + my_data_path + "/pdata_mdd/nii_xps/denoised")

ni.save(niout,  my_home_path + my_data_path + "/pdata_mdd/nii_xps/denoised/data.nii.gz")
nii.header["dim"][0] = 3
nii.header["dim"][4] = 1
