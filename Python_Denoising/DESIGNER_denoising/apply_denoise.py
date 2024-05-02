#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 10:50:01 2022

@author: omar

"""

import nibabel as ni 
import numpy as np
runfile('C:/Users/User/Mon_Drive/Python/DESIGNER_denoising/denoise.py')

nii = ni.load('C:/Users/User/Mon_Drive/Matlab/ProcessingPV360/data/20230220_Rat_Brain_3D/17/pdata_mdd/nii_xps/orig_beforeTopUp/dataBlipUp.nii.gz')

dwi = np.array(nii.dataobj)

Signal, Sigma,Nparameters = denoise(dwi, kernel=[7,7,7],patchsize=125,shrinkage='frobenius')
Sigma[np.isnan(Sigma)] = 0


niout = ni.Nifti1Image(Signal, nii.affine, nii.header)
niout.header.set_data_dtype(np.float32)
ni.save(niout, 'Mag_dn.nii')
nii.header["dim"][0] = 3
nii.header["dim"][4] = 1
niout = ni.Nifti1Image(Sigma, nii.affine, nii.header)
niout.header.set_data_dtype(np.float32)
ni.save(niout, 'Mag_sigma.nii')
niout = ni.Nifti1Image(Nparameters, nii.affine, nii.header)
niout.header.set_data_dtype(np.float32)
ni.save(niout, 'Mag_Npars.nii')



