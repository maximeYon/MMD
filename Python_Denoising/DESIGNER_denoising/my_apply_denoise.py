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

#my_home_path = "C:/Users/Administrateur/Mon_Drive/Matlab/ProcessingPV360/data/"
my_home_path = "C:/Users/User/Mon_Drive/Matlab/ProcessingPV360/data/"
my_data_path =  "EPI_RARE_MSME/15/pdata_mdd"

if not os.path.exists(my_home_path + my_data_path + "/nii_xps/orig"):
    os.makedirs(my_home_path + my_data_path + "/nii_xps/orig")
    shutil.move(my_home_path + my_data_path + "/nii_xps/data.nii.gz", my_home_path + my_data_path + "/nii_xps/orig/data.nii.gz")

nii = ni.load(my_home_path + my_data_path + "/nii_xps/orig/data.nii.gz")
dwi = np.array(nii.dataobj)

# for 2D
Signal, Sigma,Nparameters = mpden.denoise(dwi, kernel=[15,15,1],patchsize=125,shrinkage='frobenius')
# for 3D
#Signal, Sigma,Nparameters = mpden.denoise(dwi, kernel=[7,7,7],patchsize=125,shrinkage='frobenius')

niout = ni.Nifti1Image(Signal, nii.affine, nii.header)
niout.header.set_data_dtype(np.float32)

if not os.path.exists(my_home_path + my_data_path + "/nii_xps/denoised"):
    os.makedirs(my_home_path + my_data_path + "/nii_xps/denoised")

ni.save(niout,  my_home_path + my_data_path + "/nii_xps/denoised/data.nii")
nii.header["dim"][0] = 3
nii.header["dim"][4] = 1



"""  optional
Sigma[np.isnan(Sigma)] = 0
niout = ni.Nifti1Image(Sigma, nii.affine, nii.header)
niout.header.set_data_dtype(np.float32)
ni.save(niout, 'Mag_sigma.nii')
niout = ni.Nifti1Image(Nparameters, nii.affine, nii.header)
niout.header.set_data_dtype(np.float32)
ni.save(niout, 'Mag_Npars.nii')
"""
