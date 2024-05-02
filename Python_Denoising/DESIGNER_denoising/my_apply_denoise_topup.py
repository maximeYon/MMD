#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 10:50:01 2022

@author: Maxime Yon

"""

import nibabel as ni 
import numpy as np
import os
import shutil
import denoise as mpden

my_home_path = "C:/Users/Administrateur/Mon_Drive/Matlab/ProcessingPV360/data/"
# my_home_path = "C:/Users/User/Mon_Drive/Matlab/ProcessingPV360/data/"
my_data_path =  "human_sample_fresh1/34/pdata_mdd"

if not os.path.exists(my_home_path + my_data_path + "/nii_xps/orig"):
    os.makedirs(my_home_path + my_data_path + "/nii_xps/orig")
    shutil.move(my_home_path + my_data_path + "/nii_xps/data.nii.gz", my_home_path + my_data_path + "/nii_xps/orig/data.nii.gz")

nii = ni.load(my_home_path + my_data_path + "/nii_xps/orig/data.nii.gz")
dwi = np.array(nii.dataobj)
nrep = dwi.shape[3]; "dimensions start from 0"

dwiup = dwi[:,:,:,0:nrep:2];
dwidown = dwi[:,:,:,1:nrep:2];

Signalup, Sigmaup,Nparametersup = mpden.denoise(dwiup, kernel=[7,7,7],patchsize=125,shrinkage='frobenius')
Signaldown, Sigmadown,Nparametersdown = mpden.denoise(dwidown, kernel=[7,7,7],patchsize=125,shrinkage='frobenius')


nioutup = ni.Nifti1Image(Signalup, nii.affine, nii.header)
nioutup.header.set_data_dtype(np.float32)

nioutdown = ni.Nifti1Image(Signaldown, nii.affine, nii.header)
nioutdown.header.set_data_dtype(np.float32)

if not os.path.exists(my_home_path + my_data_path + "/nii_xps/denoised"):
    os.makedirs(my_home_path + my_data_path + "/nii_xps/denoised")

ni.save(nioutup,  my_home_path + my_data_path + "/nii_xps/denoised/dataUp.nii.gz")
nii.header["dim"][0] = 3
nii.header["dim"][4] = 1

ni.save(nioutdown,  my_home_path + my_data_path + "/nii_xps/denoised/dataDown.nii.gz")
nii.header["dim"][0] = 3
nii.header["dim"][4] = 1


