Pulse programs and Matlab files for acquiring and processing multidimensional dMRI on Bruker Avance spectrometers.

Written by Daniel Topgaard 20160427.
daniel.topgaard@fkem1.lu.se

Assumes familiarity with TopSpin and Matlab. 

Acquisition
 2D single-shot imaging with axisymmetric diffusion encoding.
 See fig 1 in Topgaard, Phys. Chem. Chem. Phys. 18, 8545 (2016),
 http://dx.doi.org/10.1039/c5cp07251d.
 RARE image read-out.
 Hennig et al, Magn Reson Med. 3, 823 (1986),
 http://dx.doi.org/10.1002/mrm.1910030602.

Processing
1)	Conventional diffusion tensors;
	fractional anisotropy;
	Westin's shape indices.

 2)	Isotropic and anisotropic variance of the diffusion tensor distribution;
	orientational order parameters;
 	microscopic diffusion anisotropy.
	See Lasic et al, Front. Phys. 2, 11 (2014), 
	http://dx.doi.org/10.3389/fphy.2014.00011.

 3)	Shape of the microscopic diffusion tensor (prolate, sphere, oblate).
	See, Eriksson et al., J. Chem. Phys. 142, 104201 (2015),
	http://dx.doi.org/10.1063/1.4913502.

 4)	Saupe order tensors.
	See Topgaard, Phys. Chem. Chem. Phys. 18, 8545 (2016),
	http://dx.doi.org/10.1039/c5cp07251d.

 5)	Size-shape diffusion tensor distributions.
	See de Almeida Martins and Topgaard, Phys. Rev. Lett. 116, 087601 (2016).
	http://dx.doi.org/10.1103/PhysRevLett.116.087601.
   
 6)	Size-shape-orientation diffusion tensor distributions.
	See Topgaard. J. Magn. Reson. 275, 98 (2017).
	http://dx.doi.org/10.1016/j.jmr.2016.12.007


Two versions available:
1) Avance II TopSpin 2.1.
Tested on Bruker 500 MHz with MIC-5 probe at Physical Chemistry, Lund University.

2) Avance III HD TopSpin 3.2.
Tested on Bruker 600 MHz with MIC-5 probe at Swedish NMR Center, Gothenburg.


Procedure
0) Set up Matlab paths by executing setup_paths.m.

1) Copy pulse program DT_axderare2d from
../md-dmri/acq/bruker/Avance<version>/pulseprograms
to 
/opt/topspin<version>/exp/stan/nmr/lists/pp/users.

2) Generate gradient waveforms and acquisition protocol with
bruker_axde_waveform.m and
bruker_axde_protocol.m.

3) Copy all gradient shape files g* from
../md-dmri/acq/bruker/Avance<version>/protocol
to
/opt/topspin<version>/exp/stan/nmr/lists/gp/users and
/opt/data/<user>/nmr/<dataset>/<expno>.

4) Set acquisition parameters according to the instructions in the pulse program or by copying the example data sets at
https://github.com/daniel-topgaard/md-dmri-data

5) zg (zero go: perform an acquisition)

6) Copy step1_recon_data.m and step2_run_analysis.m to
/opt/data/<user>/nmr/<dataset>/<expno>.

7) Execute step1_recon_data.m for image recon.

8) Execute step2_run_analysis.m for data processing.

9) Parameter maps in nifti and pdf format can be found in
/opt/data/<user>/nmr/<dataset>/<expno>/NII_RES/maps.

10) View the data with the GUI included in the md-dmri framework.
