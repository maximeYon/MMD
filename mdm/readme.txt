Multidimensional data management (MDM). This package contains functions
for building the local data structure.

We work with a structure often referred to as 's'. This structure has the 
following fields

s.nii_fn  - full path to a nifti file
s.xps     - the eXperimental Parameter Structure
s.mask_fn - (optional) path to a nifti file with a 3D mask
               that determined which parts of the data that contains
               actual data and not just background 

EXPERIMENTAL PARAMETER STRUCTURE (xps)

The xps has the following fields. All are in SI units. Not all fields must
be present. The model code checks that the necessary fields are present.

A) General parameters

- n:          number of images/signal values (fourth dimension). All 
              parameters need to have n entries.

- t_ex:       cumulative time from start of the experiment to the time of 
              the excitation pulse

- t_acq:      time from excitation to readout of echo, i.e., k = 0 
              (typically equals te in a SE sequence)

- intention:  a string that describes the intention of the experiment, 
              for example
                - PGSE/DTI
                - DDE/uFA
                - Angular-DDE
                - FEXI

- slice_order: not yet defined 



B) Parameters describing the total diffusion encoding effects

- b:          total b-value
- bt:         b-matrix, n x 6 (see dtd_* for format)


i) Parameters derived automatically:

- bt2:        outer product of b-matrix (fourth order tensor)
- u:          symmetry axis of bt, n x 3


ii) Fields needed where we have not yet decided on a format

- gradient waveform
- nex? number of averages per acquisition, or perhaps also noise sigma?




C) Diffusion-related parameters in a multiple diffusion encoding setting


i) The following parameters are valid for SDE, DDE et c

- mdm_delta1, mdm_delta2:            Diffusion encoding time 
                                     For SDE, we have only mdm_delta1

- mdm_capital_delta1, delta2, et c:  Time between leading edged of encoding
                                     gradients. For SDE, only 
                                     mdm_capital_delta1 is defined

- mdm_ramp_time                      Assumed to be the same throughout

- mdm_g1, mdm_g2:                    n x 3 vector, gradient amplitude

Derived parameters:

- mdm_q1, mdm_q2, ...:               q-vectors
- mdm_td1, mdm_td2, ...:             Diffusion times


ii) Parameters valid only for DDE, TDE or more 

- mde_tm12, mde_tm23, ...: diffusion-related mixing times between 
                           diffusion-encoding block 1 and 2, block 2 and 3, 
                           in a DDE, TDE, et c setting. In DDE-based FEXI, 
                           we have only mde_tm1. 

- mde_b1, mde_b2, ...:     b-value per block in a DDE, TDE et c sequence
- mde_bt1, mde_bt2, ...:   b-tensor per block in a DDE, TDE et c sequence





D) Relaxation parameters dealing with total relaxation weighting 

- te: echo time, total time with T2 weighting from excitation to readout

- tm: mixing time, total time with T1 relaxation from excitation to 
      readout

- ts: saturation recovery time, time between magnetization was zero and 
      excitation

- ti: inversion time, time between inversion of magnetization and 
      excitation

- tr: repetition time, time between excitations

In addition, we see the need for managing more complex sequence with
multiple RF pulses. We suggest to store such timings in fields called

- ste_tm: a vector of mixing times in a STE experiment

- ste_te: a vector of echo times in a STE experiment

