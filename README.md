#Multidimensional diffusion MRI 

Multidimensional diffusion MRI (MD-dMRI) is a family of conceptually related methods relying on advanced gradient modulation schemes and data processing approaches to simultaneously quantify several microstructural and dynamical properties of tissue by separating their effects on the detected MRI signal into multiple acquisition and analysis dimensions.

This GitHub repository mainly contains MATLAB software for analysis of MD-dMRI data, and also some auxiliary routines for setup of acquisition protocols, motion correction, and image visualization.

Three families of methods are currently implemented within the framework of MD-dMRI:
* Diffusion tensor distributions
* Diffusional exchange
* Diffusion and incoherent flow

The methods are described briefly [here](methods/README.md), in some more detail in two review articles,<sup>1,2</sup> and more exhaustively in the original publications cited for each method.

## How to start

Run setup_paths in the root folder to put files in the Matlab path.

To analyze data, you can use the method_name_4d_data2fit and method_name_4d_fit2param, where 'method_name' refers to one of the mehods implemented in the methods folder. Alternatively, you can use the `mdm_fit` command. See `mdm_fit --help` for an overview of how to use the command. An example is provided below: 

```matlab
mdm_fit --data input/data.nii.gz ...
    --mask tmp/mask.nii.gz ...
    --method dtd_codivide ...
    --out tmp/ ...
    --xps input/data_xps.mat;
```

The `--data` argument should point to a nifti file. The `--mask` argument is optional, but supplying e.g a brain mask may speed up the analysis by avoiding analysis of background voxels. For an overview of potential `--method` arguments, see the methods folder. The `--out` argument takes a folder/prefix as input. Finally, you need to supply data on how the experiment was performed, by giving either the `--xps` flag with the argument pointing to a .mat file holding an 'eXperiment Parameter Structure' (xps). Alternatively, for diffusion tensor distribution (dtd) methods, you can supply a b-tensor file as an argument to the `--btens` flag. 

We assume that the data has already been preprocessed before this function is called. Preprocessing should include motion and eddy current correction and potentially smoothing. These steps are also supported by the framework (see [below](#motion-and-eddy-current-correction)).

Construction of the xps or the b-tensor input can be challenging. Some help is provided [below](#construction-of-the-xps).

## GUI
The function `mgui` starts a graphical user interface that help reviewing the data quality. The left panel shows a folder structure. When opening a data file named e.g. data.nii.gz and the GUI finds an .mat file named as data_xps.mat, it will automatically load it. Then you can draw an ROI an select different methods for fitting the signal data from the dropdown menu in the 'analysis' panel to the right. 

![Image](http://markus-nilsson.github.io/md-dmri/mgui.png)

## Motion and eddy current correction

Conventional motion correction where all data is registered to a volume acquired with b = 0 s/mm<sup>2</sup> is implemented in `mdm_mec_b0`. This type of correction is not appropriate for high b-value data, where other approached such as extrapolation-based registration should be applied (see [this paper](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0141825)). This is supported in the `mdm_mec_eb` method. For an example, see below

```matlab
% Options
opt = mdm_opt;
 
% Connect to the data
s.nii_fn = fullfile(pwd, 'data.nii');
s.xps = mdm_xps_from_gdir('data_gdir.txt');
 
% Define which data to include in the reference 
% (only b-values below 1000 s/mm2)
s_ref = mdm_s_subsample(s, s.xps.b <= 1e9, fileparts(s.nii_fn), opt); 
 
% Write the elastix parameter file
p = elastix_p_affine;
p_fn = elastix_p_write(p, 'p.txt');
 
% First run a conventional coregistratin of the reference
s_ref = mdm_mec_b0(s_ref, p_fn, fileparts(s_ref.nii_fn), opt);
 
% Run an extrapolation-based registration
mdm_mec_eb(s, s_ref, p_fn);
```

## Construction of the xps
The xps hold information on how the experiment was performed. All variables are given in SI-units, and have predetermined names according to [this specification](mdm/readme.txt). Please see functions `mdm_xps_*` for help on how to construct an xps. For an example, see below

```matlab
xps = mdm_xps_from_bval_bvec('data.bval', 'data.bvec');
mdm_xps_save(xps, 'data_xps.mat')
```

If the name of the .mat-file holding the xps is related to the nifti filename as specified in `mdm_xps_fn_from_nii_fn`, the framework will automatically find the xps. 


Get acquainted with by xps structure by reading mdm/readme.txt.An extensive description of the code structure is found at http://markus-nilsson.github.io/md-dmri/.  



#References
1. D. Topgaard. Multidimensional diffusion MRI. J. Magn. Reson.,  (2017).
2. D. Topgaard. NMR methods for studying microscopic diffusion anisotropy. In: R. Valiullin (Ed.) Diffusion NMR in confined systems: Fluid transport in porous solids and heterogeneous materials, New Developments in NMR 9, Royal Society of Chemistry, Cambridge, UK (2017).

