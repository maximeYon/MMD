#!/usr/bin/env python3

import numpy as np
from dipy.data import get_fnames
from dipy.io.image import load_nifti, save_nifti
import matplotlib.pyplot as plt

from dipy.denoise.patch2self import patch2self


pathfolder = '/mnt/c/Users/User/Mon_Drive/Matlab/ProcessingPV360/data/'
pathfile = '20220524_mouse_brain_fin2/69/pdata_mdd/nii_xps/orig/data.nii.gz'


hardi_fname, hardi_bval_fname, hardi_bvec_fname = get_fnames('stanford_hardi')
data, affine = load_nifti(hardi_fname)
bvals = np.loadtxt(hardi_bval_fname)
denoised_arr = patch2self(data, bvals)