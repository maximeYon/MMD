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

To analyze data, you can use the `method_name_4d_data2fit` and `method_name_4d_fit2param` functions, where `method_name` refers to one of the mehods implemented [here](methods). Alternatively, you can use the `mdm_fit` command. See `mdm_fit --help` for an overview of how to use the command. An example is provided below: 

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
The xps holds all relevant experimental information, i.e., how the experiment was performed. All variables are given in SI-units, and have predetermined names according to [this specification](mdm/readme.txt). Please see functions `mdm_xps_*` for help on how to construct an xps. Here we give two examples. First, in the case where a single nii file is analyzed. Secondly, we show an example of how multiple nii files are merged and analyzed.

**Creating an xps based on a single nii file**

If a single nii file is the basis for the analysis, a single path is used to define its location `nii_fn = 'nii_filename.nii';`. Each nii file is usually associated with files that contain information about the diffusion encoding directions and b-values that were used. These can be used to automatically generate appropriate xps structures. To do this, call

```matlab 
xps = mdm_xps_from_gdir('gdir_filename.txt', b_delta);
``` 
or 

```matlab
xps = mdm_xps_from_bval_bvec('bval_filename.bval', 'bvec_filename.bvec', b_delta);
```

depending on what format is used. Note that the files currently do not contain information about the shape of the b-tensor, thus `b_delta` must be specified by the user. Furthermore, note that the path to gdir.txt and .bval/.bvec files can be generated based on `nii_fn` by calling either

```matlab 
gdir_fn = mdm_fn_nii2gdir(nii_fn);
```

 or 

```matlab
[bval_fn, bvec_fn] = mdm_fn_nii2bvalbvec(nii_fn);
```

The resulting xps can be saved by using `mdm_xps_save`. If the name of the .mat-file holding the xps is related to the nii filename as specified in `mdm_xps_fn_from_nii_fn`, the framework will automatically find the xps. 

**Example: xps from single nii file**

This example code requires only that the user defines the filename of the nii to create and save the corresponding xps. Note that this assumes that the gdir.txt or .bval/.bvec files are in the same folder an have standardized names.

```matlab
% Define path to nii file 
s.nii_ fn = 'path_to_nii.nii';

% We assume that the encoding is linear
b_delta   = 1;

% Get path to corresponding gdir, and generate the xps
gdir_fn   = mdm_fn_nii2gdir(s.nii_fn, b_delta);
s.xps     = mdm_xps_from_gdir(gdir_fn, b_delta);

% Or, get path to corresponding bval/bvec, and generate the xps
[bval_fn, bvec_fn] = mdm_fn_nii2bvalbvec(s.nii_fn);
s.xps     = mdm_xps_from_bval_bvec(bval_fn, bvec_fn, b_delta);

% Finally, save the xps.
xps_fn    = mdm_xps_fn_from_nii_fn(s.nii_fn);
mdm_xps_save(s.xps, xps_fn)
```


**Creating an xps based on multiple nii files**

In some cases all the necessary data cannot be acquired in a single series, which results in multiple nii files that need to be combined during the analysis. For example, encoding with linear and spherical b-tensors is currently performed in two separate series, but we wish to analyze them simultaneously. In the case when the analysis relies on multiple nii files we first create partial xps structures (according to the method for single nii files), and then we merge the nii files by using `mdm_nii_merge`, and the xps structures using `mdm_xps_merge`.

**Example: xps from multiple nii files**

This example uses alla the tools form the previous example. It assumes that full filenames are specified to all necessary nii files, and that corersponding gdir.txt or .bval/.bvec files are in the same folder and that they have standardized names (according to `mdm_fn_nii2gdir` and `mdm_fn_nii2bvalbvec`). The partial xps structures are first stored in a cell array of structures and then merged, along with the nii files. Note that this example includes only the case where gdir.txt files are used to create the xps, but can be modified to use .bval/.bvec.

```matlab
% Cell array of all files that are going to be merged
file_list = {'nii_filename_1.nii', 'nii_filename_2.nii', 'nii_filename_3.nii'};

% Corresponding b-tensor shapes. In this case: spherical, linear, and planar.
b_deltas  = [0 1 -.5];

% Define a name for the merged nii
merged_nii_fn = 'merged_nii_filename.nii';

% Loop over nii files, create partial xps structures, and store them in a cell array
for i = 1:numel(file_list)
	gdir_fn = mdm_fn_nii2gdir(file_list{i});
    xps_array{i} = mdm_xps_from_gdir(gdir_fn, [], b_deltas(i));
end

% Merge the nii files, for future use. This function also saves the nii to merged_nii_fn.
s.nii_fn = mdm_nii_merge(file_list, merged_nii_fn);

% Merge the xps structures into one. Then save the xps for future use.
s.xps = mdm_xps_merge(xps_array);
mdm_xps_save(s.xps, mdm_xps_fn_from_nii_fn(s.nii_fn));

```

Get acquainted with by xps structure by reading mdm/readme.txt.An extensive description of the code structure is found at http://markus-nilsson.github.io/md-dmri/.  


#References
1. D. Topgaard. Multidimensional diffusion MRI. J. Magn. Reson. 275, 98-113 (2017). [link](http://dx.doi.org/10.1016/j.jmr.2016.12.007)
2. D. Topgaard. NMR methods for studying microscopic diffusion anisotropy. In: R. Valiullin (Ed.) Diffusion NMR in confined systems: Fluid transport in porous solids and heterogeneous materials, New Developments in NMR 9, Royal Society of Chemistry, Cambridge, UK (2017). [link](http://dx.doi.org/10.1039/9781782623779-00226)

