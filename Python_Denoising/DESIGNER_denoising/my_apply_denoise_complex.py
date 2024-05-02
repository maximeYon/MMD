#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 10:50:01 2022

@author:Omar & Maxime

"""

import nibabel as ni 
import numpy as np
import os
import shutil
import denoise as mpden


my_home_path = "C:/Users/User/Mon_Drive/Matlab/ProcessingPV360/data/"
my_data_path =  "20230220_Rat_Brain_3D/18/pdata_mdd"

if not os.path.exists(my_home_path + my_data_path + "/nii_xps/orig"):
    os.makedirs(my_home_path + my_data_path + "/nii_xps/orig")
    shutil.move(my_home_path + my_data_path + "/nii_xps/data_real.nii.gz", my_home_path + my_data_path + "/nii_xps/orig/data_real.nii.gz")
    shutil.move(my_home_path + my_data_path + "/nii_xps/data_imag.nii.gz", my_home_path + my_data_path + "/nii_xps/orig/data_imag.nii.gz")

nii_real = ni.load(my_home_path + my_data_path + "/nii_xps/orig/data_real.nii.gz")
dwi_real = np.array(nii_real.dataobj)
SignalReal, Sigma,Nparameters = mpden.denoise(dwi_real , kernel=[7,7,7],patchsize=125,shrinkage='frobenius')

nii_imag = ni.load(my_home_path + my_data_path + "/nii_xps/orig/data_imag.nii.gz")
dwi_imag = np.array(nii_imag.dataobj)
SignalImag, Sigma,Nparameters = mpden.denoise(dwi_imag, kernel=[7,7,7],patchsize=125,shrinkage='frobenius')

SignalAbs = np.sqrt(np.square(SignalReal) + np.square(SignalImag))
nioutAbs= ni.Nifti1Image(SignalAbs, nii_real.affine, nii_real.header)
nioutAbs.header.set_data_dtype(np.float32)
if not os.path.exists(my_home_path + my_data_path + "/nii_xps/denoised"):
    os.makedirs(my_home_path + my_data_path + "/nii_xps/denoised")
    
ni.save(nioutAbs, my_home_path + my_data_path + "/nii_xps/denoised/data.nii")




