# Monte Carlo inversion of b(omega)-TR-TE encoded MRI data into
# nonparameteric D(omega)-R1-R2 distributions
#
# Code to be shared with the paper
# Massively multidimensional relaxation-diffusion correlation MRI
# O Narvaez, L Svenningsson, M Yon, A Sierra, and D Topgaard
# Front Phys special issue https://www.frontiersin.org/research-topics/21291/capturing-biological-complexity-and-heterogeneity-using-multidimensional-mri
#
# Adapted from https://github.com/daniel-topgaard/md-dmri
# Tested for Python 3.7.3

import os
import time
import math
from pathlib import Path
import numpy as np #https://numpy.org
import scipy.io as sio #https://www.scipy.org
from scipy.optimize import nnls
import nibabel as nib #https://nipy.org
import matplotlib.pyplot as plt
from pprint import pprint
import copy
import multiprocessing as mp

# Inversion options
do_multicpu = 1
range_bs = range(0,97)
range_bs = np.flipud(range_bs)

method_name = 'dtor1r2d'
dist_varnams = ('w','dpar','dperp','theta','phi','d0','rpar','rperp','r1','r2')

opt = {
    'method': method_name,
    'dmin': 5e-12,
    'dmax': 5e-9,
    'rmin': 0.1,
    'rmax': 1e5,
    'r1min': 0.2,
    'r1max': 20,
    'r2min': 2,
    'r2max': 200,
    'n_in': 200,
    'n_out': 10,
    'n_prol': 20,
    'n_darwin': 20,
    'dfuzz': .1,
    'ofuzz': .1*2*np.pi,
    'rfuzz': .1,
    'r1fuzz': .1,
    'r2fuzz': .1} # Saved later for book keeping

# Define paths
script_path = os.path.dirname(os.path.realpath(__file__))
data_path = os.path.join(script_path, 'indata')
data_name = 'data'

data_fn = os.path.join(data_path, data_name+'.nii.gz')
#mask_fn = os.path.join(data_path, data_name+'_mask_small.nii.gz')
mask_fn = os.path.join(data_path, data_name+'_mask.nii.gz')
xps_fn = os.path.join(data_path, data_name+'_xps.mat')

# Load image, mask, and experimental parameters
img = nib.load(data_fn)
hdr = img.header
ima = img.get_fdata() #array 

img = nib.load(mask_fn)
mask = np.array(img.get_fdata(), dtype=bool)
if len(mask.shape) == 2:
    mask = mask.reshape(mask.shape[0],mask.shape[1],1)
imv = ima[mask,:] #vector

# Extract relevant experimental parameters
xps_dict = sio.loadmat(xps_fn)
xps_struct = xps_dict['xps']
btomega = xps_struct[0,0]['btomega']
domega = xps_struct[0,0]['domega']
domegav = np.tile(domega[0], (1,btomega.shape[1]//6))
domegav[0,0] = 0
omega = np.cumsum(domegav,axis=1)

xps = {
    'btomega': btomega,
    'omega'  : omega,
    'tr'     : xps_struct[0,0]['tr'],
    'te'     : xps_struct[0,0]['te']}


# Define functions
# Functions to generate and manipulate DTD components
def notfin2zero(x):
    x[np.isinf(x)] = 0
    x[np.isneginf(x)] = 0
    x[np.isnan(x)] = 0

    return x

def dist_zeros(n_comp):
    # Distribution with all zeros
    zeros = np.zeros((n_comp,1))
    dist = {}
    for varnam in dist_varnams:
        dist.update( {varnam:zeros} )

    return dist

def dist_rand():
    # Distribution with random parameters within the min/max boundaries 
    dist = {
    'w'    : np.zeros((opt['n_in'],1)),
    'dpar' : opt['dmin']*np.power(opt['dmax']/opt['dmin'],np.random.rand(opt['n_in'],1)),
    'dperp': opt['dmin']*np.power(opt['dmax']/opt['dmin'],np.random.rand(opt['n_in'],1)),
    'theta': np.arccos(2*np.random.rand(opt['n_in'],1)-1),
    'phi'  : 2*np.pi*np.random.rand(opt['n_in'],1),
    'd0'   : opt['dmin']*np.power(opt['dmax']/opt['dmin'],np.random.rand(opt['n_in'],1)),
    'rpar' : opt['rmin']*np.power(opt['rmax']/opt['rmin'],np.random.rand(opt['n_in'],1)),
    'rperp': opt['rmin']*np.power(opt['rmax']/opt['rmin'],np.random.rand(opt['n_in'],1)),
    'r1'   : opt['r1min']*np.power(opt['r1max']/opt['r1min'],np.random.rand(opt['n_in'],1)),
    'r2'   : opt['r2min']*np.power(opt['r2max']/opt['r2min'],np.random.rand(opt['n_in'],1))}
 
    dpar = copy.deepcopy(dist['dpar'])
    dperp = copy.deepcopy(dist['dperp'])
    d0 = copy.deepcopy(dist['d0'])
    
    ind = d0<dpar
    d0[ind] = dpar[ind]
    ind = d0<dperp
    d0[ind] = dperp[ind]    

    for varnam in ('dpar','dperp','d0'):
        dist.update( {varnam:eval(varnam)} )
           
    return dist

def dist_sort(dist):
    # Sort distribution components in descending weight
    ind_sort = np.argsort(dist['w'])
    ind_sort = np.flipud(ind_sort)

    for varnam in dist_varnams:
        val = copy.deepcopy(dist[varnam])
        val_sort = val[ind_sort]
        dist.update( {varnam:val_sort.reshape(val_sort.size,1)} )

    return dist

def dist_extractfirst(dist,n_out):
    # Extract first (typically highest weight) components
    ind_first = np.arange(0,n_out,1)

    for varnam in dist_varnams:
        val = copy.deepcopy(dist[varnam])
        val_first = val[ind_first]
        dist.update( {varnam:val_first.reshape(val_first.size,1)} )

    return dist
    
# Functions for fitting
def xpsdist2dtomega(xps,dist):
    # Convert experimental and distribution parameters to D(omega) array 
    Nomega = xps['omega'].shape[1]
    Ncomp = dist['dpar'].shape[0]

    omega_tile = np.tile(omega, (Ncomp,1))
    dist_tile = {}
    for varnam in dist_varnams:
        val = copy.deepcopy(dist[varnam])       
        dist_tile.update( {varnam:np.tile(val.reshape(Ncomp,1), (1,Nomega))} )    
   
    dparo  = dist_tile['d0'] - (dist_tile['d0'] - dist_tile['dpar' ])/(1 + np.power(omega_tile,2)/np.power(dist_tile['rpar' ],2))
    dperpo = dist_tile['d0'] - (dist_tile['d0'] - dist_tile['dperp'])/(1 + np.power(omega_tile,2)/np.power(dist_tile['rperp'],2))
    dparo  = notfin2zero(dparo)
    dperpo = notfin2zero(dperpo)
    dperpo[dperpo==0] = opt['dmin']
    dparo[dparo==0]   = opt['dmin']
    
    diso = (dparo + 2*dperpo)/3
    ddelta = notfin2zero((dparo-dperpo)/(3*diso))

    xcos = np.cos(dist_tile['phi'])*np.sin(dist_tile['theta'])
    ycos = np.sin(dist_tile['phi'])*np.sin(dist_tile['theta'])
    zcos = np.cos(dist_tile['theta'])

    dxx = diso*(1 + ddelta*(3*xcos*xcos - 1))
    dxy = diso*(0 + ddelta*(3*xcos*ycos - 0))
    dxz = diso*(0 + ddelta*(3*xcos*zcos - 0))
    dyy = diso*(1 + ddelta*(3*ycos*ycos - 1))
    dyz = diso*(0 + ddelta*(3*ycos*zcos - 0))
    dzz = diso*(1 + ddelta*(3*zcos*zcos - 1))

    sqrt2 = np.sqrt(2)
    dtomega = np.concatenate((dxx,dyy,dzz,sqrt2*dxy,sqrt2*dxz,sqrt2*dyz),axis=1)

    return dtomega
    
def kernel(xps,dist):
    # Generate inversion kernel from experimental and distribution parameters
    dtomega = xpsdist2dtomega(xps,dist)
    Kd = np.exp(-np.matmul(xps['btomega'], dtomega.T))
    Kr1 = 1-np.exp(-np.matmul(xps['tr'], dist['r1'].T))
    Kr2 = np.exp(-np.matmul(xps['te'], dist['r2'].T))
    K = Kd*Kr1*Kr2
    return K

def s2w_inversion(s,xps,dist):
    # Inversion from signal to component weights. Outputs sorted distribution
    smax = np.max(s)
    dist.update( {'w': smax*fnnls_Ks(kernel(xps,dist),s/smax)})
    dist = dist_sort(dist)     
    return dist

def fnnls_Ks(K,s):
    # Preparation for convenient switching between inversion algorithms

    # Fast nnls for large matrices
    # KtK = K.T.dot(K)
    # Kts = K.T.dot(s)
    # w = fnnls(KtK, Kts)

    # Standard nnls
    w, rnorm = nnls(K, s) 

    return w

def signal2dist(s):
    # Monte Carlo inversion of signal to sorted distribution ready for saving
    if np.max(s)>0:
        dist = dist_zeros(opt['n_in'])

        # Proliferation
        for count in range(opt['n_prol']):
            ind_zero = dist['w'] == 0
            dist_new = dist_rand()
            for varnam in dist_varnams:
                val = copy.deepcopy(dist[varnam])
                val = val.reshape(val.size,1)
                val_new = copy.deepcopy(dist_new[varnam])
                val[ind_zero] = val_new[ind_zero]
                dist.update( {varnam:val} )
            
            dist = s2w_inversion(s,xps,dist)
           
        # Darwinian mutation and extinction
        for count in range(opt['n_darwin']): 
            ind_nonzero = np.argwhere(dist['w']>0)
            n_max = np.min([opt['n_out'], np.size(ind_nonzero)])

            # Assure that high-weight components produce more offspring
            ind_keep = np.floor((n_max-1)*np.power(np.linspace(0,1,opt['n_in']-n_max),3))
            ind_keep = ind_keep.astype(int)
            ind_replace = np.arange(n_max,opt['n_in'],1)
            Nreplace = ind_replace.size
            
            # Temporarily copy distribution in dict to arrays
            w     = copy.deepcopy(dist['w'])
            dpar  = copy.deepcopy(dist['dpar'])
            dperp = copy.deepcopy(dist['dperp'])
            theta = copy.deepcopy(dist['theta'])
            phi   = copy.deepcopy(dist['phi'])
            d0    = copy.deepcopy(dist['d0'])
            rpar  = copy.deepcopy(dist['rpar'])
            rperp = copy.deepcopy(dist['rperp'])
            r1    = copy.deepcopy(dist['r1'])
            r2    = copy.deepcopy(dist['r2'])

            dpar_new  = dpar
            dperp_new = dperp
            theta_new = theta
            phi_new   = phi
            d0_new    = d0
            rpar_new  = rpar
            rperp_new = rperp
            r1_new    = r1
            r2_new    = r2

            # Replacement of the weakest by mutated fittest
            dpar[ind_replace]  = dpar_new[ind_keep]*(1+opt['dfuzz']*np.random.randn(Nreplace,1))
            dperp[ind_replace] = dperp_new[ind_keep]*(1+opt['dfuzz']*np.random.randn(Nreplace,1))
            theta[ind_replace] = theta_new[ind_keep]+opt['ofuzz']*np.random.randn(Nreplace,1)
            phi[ind_replace]   = phi_new[ind_keep]+opt['ofuzz']*np.random.randn(Nreplace,1)
            d0[ind_replace]    = d0_new[ind_keep]*(1+opt['dfuzz']*np.random.randn(Nreplace,1))
            rpar[ind_replace]  = rpar_new[ind_keep]*(1+opt['rfuzz']*np.random.randn(Nreplace,1))
            rperp[ind_replace] = rperp_new[ind_keep]*(1+opt['rfuzz']*np.random.randn(Nreplace,1))
            r1[ind_replace]    = r1_new[ind_keep]*(1+opt['r1fuzz']*np.random.randn(Nreplace,1))
            r2[ind_replace]    = r2_new[ind_keep]*(1+opt['r2fuzz']*np.random.randn(Nreplace,1))

            # Clamping
            dpar[dpar>opt['dmax']] = opt['dmax']
            dpar[dpar<opt['dmin']] = opt['dmin']
            dperp[dperp>opt['dmax']] = opt['dmax']
            dperp[dperp<opt['dmin']] = opt['dmin']
            d0[d0>opt['dmax']] = opt['dmax']
            d0[d0<opt['dmin']] = opt['dmin']
            ind = d0<dpar
            d0[ind] = dpar[ind]
            ind = d0<dperp
            d0[ind] = dperp[ind]
            rpar[rpar>opt['rmax']] = opt['rmax']
            rpar[rpar<opt['rmin']] = opt['rmin']
            rperp[rperp>opt['rmax']] = opt['rmax']
            rperp[rperp<opt['rmin']] = opt['rmin']
            r1[r1>opt['r1max']] = opt['r1max']
            r1[r1<opt['r1min']] = opt['r1min']
            r2[r2>opt['r2max']] = opt['r2max']
            r2[r2<opt['r2min']] = opt['r2min']

            for varnam in dist_varnams:
                dist.update( {varnam:eval(varnam)} )

            dist = s2w_inversion(s,xps,dist)
           
        if opt['n_out']<opt['n_in']:
            # Select top fittest and re-weight
            dist = dist_extractfirst(dist,opt['n_out'])
            dist = s2w_inversion(s,xps,dist)
                   
    else:
        dist = dist_zeros(opt['n_out'])

    # Null distribution parameters for components with zero weight and
    # convert to single precision to save disk space 
    ind_zero = dist['w'] == 0
    for varnam in dist_varnams:
        val = np.float32(copy.deepcopy(dist[varnam]))
        val[ind_zero] = 0
        dist.update( {varnam:val.reshape(val.size,1)} )

    # lines = plt.plot(range(0,s.shape[0]), s, 'ro', range(0,s.shape[0]), np.matmul(kernel(xps,dist),dist['w']), 'k.')
    # plt.show() 
   
    # Copy distribution dict for saving
    fit = copy.deepcopy(dist)

    return fit


# Prepare for looping over voxels and bootstrapping
n_vox = imv.shape[0]
n_samp = imv.shape[1]
imv_full = copy.deepcopy(imv)
xps_full = copy.deepcopy(xps)
xps_varnams = ('btomega','tr','te')

# Loop over bootstraps
if __name__ == '__main__':
    for count_bs in range_bs:
        tic = time.perf_counter()

        # Bootstrap resampling
        ind_resample = np.floor(np.random.rand(n_samp)*n_samp)
        ind_resample = ind_resample.astype(int)
        # ind_resample = range(0,n_samp) #Skip resampling for debugging
        imv = imv_full[:,ind_resample]
    
        for varnam in xps_varnams:
            xps.update( {varnam:copy.deepcopy(xps_full[varnam][ind_resample,:])} )
    
        # Loop over voxels
        range_vox = range(n_vox)

        if do_multicpu:
            # Multi CPUs
            pool = mp.Pool(mp.cpu_count())
            fit_list = pool.map(signal2dist, [imv[count_vox,:] for count_vox in range_vox])
            pool.close()
        else:
            #Single CPU
            fit_list = [None] * imv.shape[0]
            for count_vox in range_vox:
                fit_list[count_vox] = signal2dist(imv[count_vox,:])

        # Save individual bootstraps to file
        (head, tail) = os.path.split(data_path)
        out_path = os.path.join(head, 'distributions', 'bootstraps', str(count_bs))
        Path(out_path).mkdir(parents=True, exist_ok=True)
        mfs_fn = os.path.join(out_path, 'mfs.npz')
        np.savez_compressed(mfs_fn,fit_list=fit_list,mask=mask,opt=opt,ind_resample=ind_resample)

        toc = time.perf_counter()
        print(f"Processed in {(toc - tic):0.4f} seconds, {imv.shape[0]/(toc - tic):0.4f} voxels/second")