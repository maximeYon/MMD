#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 24/04/2023

@author: Maxime Yon

"""

import nibabel as ni 
import numpy as np
import numpy.matlib
import os
import sys
import denoise as mpden

# my_home_path = "C:/Users/Administrateur/Mon_Drive/Matlab/ProcessingPV360/data/"
my_data_path = str(sys.argv[1]);
# my_data_path = "ON-81-mbti-pilot-3d/17"

nii = ni.load(my_data_path + "/data.nii.gz")
dwi = np.array(nii.dataobj)
nrep = dwi.shape[3]; "dimensions start from 0"

dwiup = dwi[:,:,:,0:nrep:2];
dwidown = dwi[:,:,:,1:nrep:2];

Signalup, Sigmaup,Nparametersup = mpden.denoise(dwiup, kernel=[15,15,5],patchsize=125,shrinkage='frobenius')
Signaldown, Sigmadown,Nparametersdown = mpden.denoise(dwidown, kernel=[15,15,5],patchsize=125,shrinkage='frobenius')

#Signalup, Sigmaup,Nparametersup = mpden.denoise(dwiup, kernel=[7,7,7],patchsize=125,shrinkage='frobenius')
#Signaldown, Sigmadown,Nparametersdown = mpden.denoise(dwidown, kernel=[7,7,7],patchsize=125,shrinkage='frobenius')

# Rician Noise correction corresponding to:
#run.command('mrcalc working.mif 2 -pow lowbnoisemap.mif 2 -pow -sub -abs -sqrt - | mrcalc - -finite - 0 -if dwirc.mif')
Signalup = np.power(Signalup,2);
Sigmaup = np.power(Sigmaup,2);
Sigmaup = numpy.matlib.repmat(Sigmaup,1,1,1,nrep/2)
Signalup = np.sqrt(np.abs(Signalup-Sigmaup));
Signalup = np.multiply(Signalup,np.isfinite(Signalup));

Signaldown = np.power(Signaldown,2);
Sigmadown = np.power(Sigmadown,2);
Sigmadown = numpy.matlib.repmat(Sigmadown,1,1,1,nrep/2)
Signaldown = np.sqrt(np.abs(Signaldown-Sigmadown));
Signaldown = np.multiply(Signaldown,np.isfinite(Signaldown));
# End of Rician Noise correction

nioutup = ni.Nifti1Image(Signalup, nii.affine, nii.header)
nioutup.header.set_data_dtype(np.float32)

nioutdown = ni.Nifti1Image(Signaldown, nii.affine, nii.header)
nioutdown.header.set_data_dtype(np.float32)

if not os.path.exists(my_data_path + "/denoised"):
       os.makedirs(my_data_path + "/denoised")

ni.save(nioutup,  my_data_path + "/denoised/dataUp.nii.gz")
nii.header["dim"][0] = 3
nii.header["dim"][4] = 1

ni.save(nioutdown,  my_data_path + "/denoised/dataDown.nii.gz")
nii.header["dim"][0] = 3
nii.header["dim"][4] = 1