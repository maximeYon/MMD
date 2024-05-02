# Conversion from nonparameteric D(omega)-R1-R2 distributions
# to parameter maps
#
# Code to be shared with the paper
# Massively multidimensional relaxation-diffusion correlation MRI
# O Narvaez, L Svenningsson, M Yon, A Sierra, and D Topgaard
# Front Phys special issue https://www.frontiersin.org/research-topics/21291/capturing-biological-complexity-and-heterogeneity-using-multidimensional-mri
#
# Adapted from https://github.com/daniel-topgaard/md-dmri
# Tested for Python 3.7.3

import os
import glob
import time
import math
from pathlib import Path
import numpy as np #https://numpy.org
np.seterr(divide='ignore', invalid='ignore') # turn off warning for divide by zero
import scipy.io as sio #https://www.scipy.org
from scipy.optimize import nnls
import nibabel as nib #https://nipy.org
import matplotlib.pyplot as plt #https://matplotlib.org/#
from pprint import pprint
import multiprocessing as mp

do_multicpu = 1

# Define paths
script_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(script_path, 'indata')
data_name = 'data'
method_name = 'dtor1r2d'

data_fn = os.path.join(data_path, data_name+'.nii.gz')
xps_fn = os.path.join(data_path, data_name+'_xps.mat')

# Define frequencies for evaluating the diffusion dimensions
omega = 2*np.pi*50 #For most maps
omega_high = 2*np.pi*150 #High-frequency for restriction

# Define color bar limits
clim = {'mdiso': [0, 2e-9],
    'msddelta': [0, 1],
    'mr1': [0, 10],
    'mr2': [0, 150],
    'cvdisosddelta': .27e-9*np.array([-1, 1]),
    'cvdisor1': 2e-9*np.array([-1, 1]),
    'cvdisor2': 2e-8*np.array([-1, 1]),
    'cvr1r2': 100*np.array([-1, 1]),
    'dmdisodnu': 2e-12*np.array([-1, 1]),
    'dmsddeltadnu': 1e-3*np.array([-1, 1]),
    'dvdisodnu': .5e-21*np.array([-1, 1]),
    'dvsddeltadnu': .15e-3*np.array([-1, 1]),
    'dcvdisosddeltadnu': 0.274e-12*np.array([-1, 1])}
# Define bins in diso-sddelta space
binlim = {'diso_min' : np.array([  0,  0, 1])*1e-9,
    'diso_max'       : np.array([  1,  1, 5])*1e-9,
    'sddelta_min'    : np.array([.25,  0, 0]),
    'sddelta_max'    : np.array([  1,.25, 1])}
binlim.update( {'Nbins': len(binlim['diso_min'])} )

# Find completed bootstraps
(head, tail) = os.path.split(data_path)
mfs_fn_pattern = os.path.join(head, 'distributions', 'bootstraps','*','mfs.npz')
mfs_fns=glob.glob(mfs_fn_pattern)
Nbs = len(mfs_fns)
# Nbs = 8 # To check the denoising effect from multiple bootstrap repetitions


# Functions

# From fits to stats
def mfs_fn2paramsprim(mfs_fn):
    # Load primary fit parameters
    # Rearrange from list of per-voxel dict to dict with per-parameter array (component x voxel)

    data = np.load(mfs_fn,allow_pickle=True)
    fit_list = data['fit_list']

    varnams = list(fit_list[1].keys())
    Nout = len(fit_list[1][varnams[1]]) # Number of distribution components
    Nvoxels = len(fit_list) # Number of voxels

    # Initialize fit dict as numpy arrays components x voxels
    paramsprim = {}
    for varnam in varnams:
        value = np.zeros([Nout,Nvoxels])
        paramsprim.update( {varnam:value} )

    # Loop over voxels
    counter = range(Nvoxels)
    for nvoxel in counter:
        fit = fit_list[nvoxel]
        if fit:
            for varnam in varnams:
                value = fit[varnam]
                paramsprim[varnam][:,nvoxel] = value.flatten()

    return paramsprim
    
def paramsprim2paramssec(paramsprim):
    # Convert from primary to secondary parameters

    w = paramsprim['w']
    dpar = paramsprim['dpar']
    dperp = paramsprim['dperp']
    theta = paramsprim['theta']
    phi = paramsprim['phi']
    d0 = paramsprim['d0']
    rpar = paramsprim['rpar']
    rperp = paramsprim['rperp']
    r1 = paramsprim['r1']
    r2 = paramsprim['r2']

    # Convert from Lorentzian parameters to DTD at specific omega
    dparo = d0 - (d0 - dpar)/(1 + np.power(omega,2)/np.power(rpar,2))
    dperpo = d0 - (d0 - dperp)/(1 + np.power(omega,2)/np.power(rperp,2))

    diso = (dparo + 2*dperpo)/3
    ddelta = notfin2zero((dparo-dperpo)/3/diso)
    sddelta = np.power(ddelta,2)
    xcos = np.cos(phi)*np.sin(theta)
    ycos = np.sin(phi)*np.sin(theta)
    zcos = np.cos(theta)
    dxx = diso*(1 + ddelta*(3*xcos*xcos - 1))
    dyy = diso*(1 + ddelta*(3*ycos*ycos - 1))
    dzz = diso*(1 + ddelta*(3*zcos*zcos - 1))
    dxy = diso*(0 + ddelta*(3*xcos*ycos - 0))
    dxz = diso*(0 + ddelta*(3*xcos*zcos - 0))
    dyz = diso*(0 + ddelta*(3*ycos*zcos - 0))

    # Repeat for diso and sddelta at high omega
    dparo = d0 - (d0 - dpar)/(1 + np.power(omega_high,2)/np.power(rpar,2))
    dperpo = d0 - (d0 - dperp)/(1 + np.power(omega_high,2)/np.power(rperp,2))

    diso_high = (dparo + 2*dperpo)/3
    ddelta_high = notfin2zero((dparo-dperpo)/3/diso_high)
    sddelta_high = np.power(ddelta_high,2)

    varnams = ('w','diso','sddelta','dxx','dyy','dzz','dxy','dxz','dyz','diso_high','sddelta_high','r1','r2')
    paramssec = {}
    for varnam in varnams:
        paramssec.update( {varnam:eval(varnam)} )
     
    return paramssec

def paramssec2stats(paramssec,mask_bin):
    # Convert from secondary parameters to statistical descriptors
    # for one bin

    # Synthetic D, T1, and T2 weighted images
    param = paramssec
    stats = {'s0':np.sum(mask_bin*param['w'],axis=0)}
    stats.update( {'b0500':np.sum(mask_bin*param['w']*np.exp(-.5e9*param['diso']),axis=0)} )
    stats.update( {'b1000':np.sum(mask_bin*param['w']*np.exp(-1e9*param['diso']),axis=0)} )
    stats.update( {'b2000':np.sum(mask_bin*param['w']*np.exp(-2e9*param['diso']),axis=0)} )
    stats.update( {'b3000':np.sum(mask_bin*param['w']*np.exp(-3e9*param['diso']),axis=0)} )
    stats.update( {'tr100':np.sum(mask_bin*param['w']*(1-np.exp(-100e-3*param['r1'])),axis=0)} )
    stats.update( {'tr200':np.sum(mask_bin*param['w']*(1-np.exp(-200e-3*param['r1'])),axis=0)} )
    stats.update( {'tr500':np.sum(mask_bin*param['w']*(1-np.exp(-500e-3*param['r1'])),axis=0)} )
    stats.update( {'te010':np.sum(mask_bin*param['w']*np.exp(-10e-3*param['r2']),axis=0)} )
    stats.update( {'te020':np.sum(mask_bin*param['w']*np.exp(-20e-3*param['r2']),axis=0)} )
    stats.update( {'te050':np.sum(mask_bin*param['w']*np.exp(-50e-3*param['r2']),axis=0)} )
    # Means
    varnams = ('diso','sddelta','r1','r2',
        'dxx','dyy','dzz','dxy','dxz','dyz',
        'diso_high','sddelta_high')
    for varnam in varnams:
        stats.update( {'m'+varnam:np.sum(mask_bin*param['w']*param[varnam],axis=0)/stats['s0']} )
    # Variances
    varnams = ('diso','sddelta','r1','r2','diso_high','sddelta_high')
    for varnam in varnams:
        stats.update( {'v'+varnam:np.sum(mask_bin*param['w']*np.power(param[varnam]-stats['m'+varnam],2),axis=0)/stats['s0']} )
    # Covariances
    varnams1=(   'diso','diso','diso','sddelta','sddelta','r1')
    varnams2=('sddelta',  'r1',  'r2',     'r1',     'r2','r2')
    for count in range(len(varnams1)):
        varnam1 = varnams1[count]
        varnam2 = varnams2[count]
        stats.update( {'cv'+varnam1+varnam2:np.sum(mask_bin*param['w']*(param[varnam1]-stats['m'+varnam1])*(param[varnam2]-stats['m'+varnam2]),axis=0)/stats['s0']} )
    # Frequency-dependencies
    varnams = ('diso','sddelta')
    for varnam in varnams:
        # Means
        value_high = np.sum(mask_bin*param['w']*param[varnam+'_high'],axis=0)/stats['s0']
        value = np.sum(mask_bin*param['w']*param[varnam],axis=0)/stats['s0']
        stats.update( {'dm'+varnam+'dnu':(value_high-value)/((omega_high-omega)/(2*np.pi))} )
        # Variances
        value_high = np.sum(mask_bin*param['w']*np.power(param[varnam+'_high']-stats['m'+varnam+'_high'],2),axis=0)/stats['s0']
        value = np.sum(mask_bin*param['w']*np.power(param[varnam]-stats['m'+varnam],2),axis=0)/stats['s0']
        stats.update( {'dv'+varnam+'dnu':(value_high-value)/((omega_high-omega)/(2*np.pi))} )
    # Covariance
    varnam1 = 'diso'
    varnam2 = 'sddelta'
    value_high = np.sum(mask_bin*param['w']*(param[varnam1+'_high']-stats['m'+varnam1+'_high'])*(param[varnam2+'_high']-stats['m'+varnam2+'_high']),axis=0)/stats['s0']
    value = np.sum(mask_bin*param['w']*(param[varnam1]-stats['m'+varnam1])*(param[varnam2]-stats['m'+varnam2]),axis=0)/stats['s0']
    stats.update( {'dcv'+varnam1+varnam2+'dnu':(value_high-value)/((omega_high-omega)/(2*np.pi))} )

    return stats
    
def paramssec2binstats(paramssec,binlim):
    # Convert from secondary parameters to statistical descriptors
    # as list over all bins

    binstats = [None] * (binlim['Nbins']+1)
    # Global
    binstats[0] = paramssec2stats(paramssec,1)
    # Bins
    Nout = paramssec['w'].shape[0]
    Nvoxels = paramssec['w'].shape[1]
    for nbin in range(binlim['Nbins']):
        mask_bin = np.zeros([Nout,Nvoxels,4], dtype=bool)
        mask_bin[:,:,0] = paramssec['diso']   <=binlim['diso_max'][nbin]
        mask_bin[:,:,1] = paramssec['diso']   >=binlim['diso_min'][nbin]
        mask_bin[:,:,2] = paramssec['sddelta']<=binlim['sddelta_max'][nbin]
        mask_bin[:,:,3] = paramssec['sddelta']>=binlim['sddelta_min'][nbin]
        mask_bin = np.all(mask_bin,axis=2)
        binstats[nbin+1] = paramssec2stats(paramssec,mask_bin)

    return binstats

def mfs_fn2binstats(mfs_fn):
    # From file names to list of per-bin stats

    paramsprim = mfs_fn2paramsprim(mfs_fn)
    paramssec = paramsprim2paramssec(paramsprim)
    binstats = paramssec2binstats(paramssec,binlim)

    return binstats

# Clamping
def notfin2zero(x):
    x[np.isinf(x)] = 0
    x[np.isneginf(x)] = 0
    x[np.isnan(x)] = 0
    return x

def clamp01(x):
    x[x<0] = 0
    x[x>1] = 1
    x[~np.isfinite(x)] = 0
    return x

def clamp_m1to1(x):
    x[x<-1] = -1
    x[x>1] = 1
    x[~np.isfinite(x)] = 0
    return x

# Color maps
def cind2rgb(cind):
    # Convert from color index 0<cind<1 to RGB blue-bright green-red map
    cind = clamp01(cind)

    r = -1 + 3*cind
    g =  1.5 - 3*np.abs(cind-.5)
    b =  2 - 3*cind

    r = clamp01(r)
    g = clamp01(g)
    b = clamp01(b)
    return r,g,b

def cind2rgb_hotcold(cind):
    # Convert from color index -1<cind<1 to RGB hotcold map
    cind = clamp_m1to1(cind)

    r =  0 + 8/3*cind
    g = -1 + 8/3*cind
    b = -3 + 8/2*cind

    ind_neg = cind < 0

    b[ind_neg] =  0 - 8/3*cind[ind_neg]
    g[ind_neg] = -1 - 8/3*cind[ind_neg]
    r[ind_neg] = -3 - 8/2*cind[ind_neg]

    r = clamp01(r)
    g = clamp01(g)
    b = clamp01(b)
    return r,g,b

def rgb_norm(r,g,b):
    # Normalize RGB values
    rgb_norm = np.zeros([3,r.size])
    rgb_norm[0,:] = r
    rgb_norm[1,:] = g
    rgb_norm[2,:] = b
    rgb_norm = np.amax(rgb_norm,axis=0)
    r = r/rgb_norm
    g = g/rgb_norm
    b = b/rgb_norm
    return r,g,b

# From stats to niftis
def im1d_to_im3d(im1d,mask):
    # Convert image from 1D vector to 3D array
    im3d = np.zeros([*mask.shape])
    im3d[mask,...] = im1d
    return im3d

def brightrgb_to_nii(nii_fn,bright,r,g,b,mask):
    # Convert brightness and RGB image vectors to 3D arrays and save as nifti

    r = im1d_to_im3d(r,mask)
    g = im1d_to_im3d(g,mask)
    b = im1d_to_im3d(b,mask)
    bright = im1d_to_im3d(bright,mask)

    im3d_rgb = np.zeros([*mask.shape, 3])
    im3d_rgb[:,:,:,0] = bright*r
    im3d_rgb[:,:,:,1] = bright*g
    im3d_rgb[:,:,:,2] = bright*b
    im3d_rgb = np.uint8(clamp01(im3d_rgb)*255)
    ras_pos = im3d_rgb
    # ras_pos is a 4-d numpy array, with the last dim holding RGB
    shape_3d = ras_pos.shape[0:3]
    rgb_dtype = np.dtype([('R', 'u1'), ('G', 'u1'), ('B', 'u1')])
    ras_pos = ras_pos.copy().view(dtype=rgb_dtype).reshape(shape_3d)  # copy used to force fresh internal structure
    # Hack that should be done more elegantly: First save nifti with default header to get right data format for RGB
    ni_img = nib.Nifti1Image(ras_pos, np.eye(4))
    nib.save(ni_img, nii_fn)
    # ... and then reload and save with dimensions from original header
    resave_nii(nii_fn,hdr)

def resave_nii(nii_fn,hdr):
    # ... and then reload and save with dimensions from original header
    img_new = nib.load(nii_fn)
    hdr_new = img_new.header
    data_new = np.asanyarray(img_new.dataobj)
    fieldnams = ('pixdim','vox_offset','scl_slope','scl_inter','descrip','aux_file',
                'qform_code','sform_code','quatern_b','quatern_c','quatern_d',
                'qoffset_x','qoffset_y','qoffset_z','srow_x','srow_y','srow_z',
                'intent_name','magic')
                
    for fieldnam in fieldnams:
        hdr_new[fieldnam] = hdr[fieldnam]

    ni_img = nib.Nifti1Image(data_new, affine=None, header=hdr_new)
    nib.save(ni_img, nii_fn)

# End functions


# Load fits, convert to stats list over bins, and compile in list over bootstraps
if __name__ == '__main__':
    if do_multicpu:
        # Multi CPUs
        pool = mp.Pool(mp.cpu_count()) 
        statslist_bslist = pool.map(mfs_fn2binstats, [mfs_fns[nbs] for nbs in range(Nbs)])
        pool.close()
    else:
        # Single CPU
        statslist_bslist = [None] * Nbs
        for nbs in range(Nbs):
            statslist_bslist[nbs] = mfs_fn2binstats(mfs_fns[nbs])



    data = np.load(mfs_fns[0],allow_pickle=True)
    mask = data['mask']

    Nbs = len(statslist_bslist)
    Nbins = len(statslist_bslist[0])
    varnam = 's0'
    Nvoxels = len(statslist_bslist[0][0][varnam])

    # Calculate medians over bootstraps
    binstats = [None] * Nbins
    for nbin in range(Nbins):
        binstats[nbin] = {}
        for varnam in statslist_bslist[0][nbin].keys():
            val_array = np.full([Nvoxels,Nbs],np.nan)
            for nbs in range(Nbs):
                val_array[:,nbs] = statslist_bslist[nbs][nbin][varnam]
            val = np.nanmedian(val_array,axis=1)
            binstats[nbin].update( {varnam:val} )

    # Save parameter maps as nifti
    (head, tail) = os.path.split(data_path)
    # out_path = os.path.join(head, 'pmaps_py')
    out_path = os.path.join(head,'maps','Nbs_'+str(Nbs))
    Path(out_path).mkdir(parents=True, exist_ok=True)
    pmap_ext = '.nii.gz'

    # Load original nifti to get header
    img = nib.load(data_fn)
    hdr = img.header

    # 1. Gray scale maps for positive scalar parameters
    # Convert stats to parameter map dict
    stats = binstats[0]
    shape = [*mask.shape]
    pmap_dict = {}
    pmap_prefix = 'dtor1r2d_'
    varnams = ('s0','b0500','b1000','b2000','b3000',
    'tr100','tr200','tr500',
    'te010','te020','te050',
    'mdiso','msddelta','vdiso','vsddelta',
    'mr1','mr2','vr1','vr2')

    for varnam in varnams:
        pmapnam = pmap_prefix+varnam
        pmap_dict.update( {pmapnam:np.zeros(shape)} )
        pmap_dict[pmapnam][mask,...] = np.float32(notfin2zero(stats[varnam]))

    for pmapnam in pmap_dict.keys():
        im3d = pmap_dict[pmapnam]
        nii_fn = os.path.join(out_path, pmapnam+pmap_ext)
        # img_out = nib.Nifti1Image(im3d, affine=None, header=hdr)
        img_out = nib.Nifti1Image(im3d, np.eye(4))
        nib.save(img_out, nii_fn)
        resave_nii(nii_fn,hdr)

    # 2. Hotcold maps for scalar parameters spanning negative and positive values
    varnams = ('cvdisosddelta','cvdisor1','cvdisor2','cvr1r2',
    'dmdisodnu','dmsddeltadnu','dvdisodnu','dvsddeltadnu','dcvdisosddeltadnu')

    for varnam in varnams:
        pmapnam = pmap_prefix+varnam
        nii_fn = os.path.join(out_path, pmapnam+pmap_ext)
        val = binstats[0][varnam]
        cind = 2*val/(clim[varnam][1]-clim[varnam][0])
        r,g,b = cind2rgb_hotcold(cind)
        bright = np.ones([*r.shape])
        brightrgb_to_nii(nii_fn,bright,r,g,b,mask)

    # 3. Fractions (brightness) and means (colors) for bin parameters
    for nbin in range(1,Nbins):
        bright = notfin2zero(binstats[nbin]['s0']/binstats[0]['s0'])
        varnams = ('mdiso','msddelta','mr1','mr2','dmdisodnu','dmsddeltadnu')
        for varnam in varnams:
            pmapnam = pmap_prefix+varnam+'_bin'+str(nbin)
            nii_fn = os.path.join(out_path, pmapnam+pmap_ext)

            val = binstats[nbin][varnam]
            cind = (val-clim[varnam][0])/(clim[varnam][1]-clim[varnam][0])
            r,g,b = cind2rgb(cind)

            brightrgb_to_nii(nii_fn,bright,r,g,b,mask)
        
        pmapnam = pmap_prefix+'mdii_bin'+str(nbin)
        nii_fn = os.path.join(out_path, pmapnam+pmap_ext)

        r = binstats[nbin]['mdxx']
        g = binstats[nbin]['mdyy']
        b = binstats[nbin]['mdzz']  
        r,g,b = rgb_norm(r,g,b)
        brightrgb_to_nii(nii_fn,bright,r,g,b,mask)    

    # 4. Special treatment of RGB fractions map
    pmapnam = pmap_prefix+'fractions'
    nii_fn = os.path.join(out_path, pmapnam+pmap_ext)

    r = binstats[1]['s0']
    g = binstats[2]['s0']
    b = binstats[3]['s0']  
    bright = np.ones([*r.shape])
    r,g,b = rgb_norm(r,g,b)
    brightrgb_to_nii(nii_fn,bright,r,g,b,mask)