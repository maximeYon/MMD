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

All parameters in the xps relating to different acquisitions are stored as
variables of size n x m, where n is the number of image volumes and m is 
the number of parameters. 

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

- a_ind:      Averaging index. After averaging (arithmetic or geometric),
              there will be max(a_ind) number of images left.

- s_ind:      Series index. Refers to data acquired in different series, 
              for example, with different prescans (e.g. gain adjustment).
              Can also index acquisitions with different echo times et c.
              

Potentially but not necessarily automatically calculated:

- b_ind:      Indexes measurements according to total b-values

- bd_ind:     Indexes measurements according to b-anisotropy (b_delta)

- be_ind:     Indexes measurements according to b-asymmetry (b_eta)

- br_ind:     Indexes measurements according to b-tensor rotations



B) Parameters describing the total diffusion encoding effects

- b:          total b-value
- bt:         b-matrix, n x 6 (see dtd_* for format)
- alpha:      flow compensation factor (see Ahlgren 16)
- alpha2:     flow attenuation factor (see Ahlgren 16)


i) Parameters derived automatically:

- bt2:        outer product of b-matrix (fourth order tensor)
- u:          symmetry axis of bt, n x 3

- b_delta:    anisotropy of b-tensor
- b_eta:      asymmetry of b-tensor
- b_s:        spherical component of b-tensor (Martins et al, PRL 16)
- b_p:        planar component of b-tensor (Martins et al, PRL 16)
- b_l:        linear component of b-tensor (Martins et al, PRL 16)

- b_paszz:    b-eigenvalue furthest from the mean ("symmetry axis").
- b_pasyy:    b-eigenvalue closest from the mean ("symmetry axis")
- b_pasxx:    b-eigenvalue that is not zz or yy (Eriksson et al, JCP 15)
- b_alpha:    Euler angle of the b-tensor, from LAB TO PAS
- b_beta:     Euler angle of the b-tensor
- b_gamma:    Euler angle of the b-tensor
    


ii) Fields needed where we have not yet decided on a format

- gradient waveform
- nex? number of averages per acquisition, or perhaps also noise sigma?




C) Diffusion-related parameters in a multiple diffusion encoding setting


i) The following parameters are valid for SDE, DDE et c

- mde_delta1, mde_delta2:            Diffusion encoding time 
                                     For SDE, we have only mde_delta1

- mde_capital_delta1, delta2, et c:  Time between leading edged of encoding
                                     gradients. For SDE, only 
                                     mde_capital_delta1 is defined

- mde_ramp_time                      Assumed to be the same throughout

- mde_g1, mde_g2:                    n x 3 vector, gradient amplitude

Derived parameters:

- mde_q1, mde_q2, ...:               q-vectors
- mde_td1, mde_td2, ...:             Diffusion times


ii) Parameters valid only for DDE, TDE or more 

- mde_tm12, mde_tm23, ...: diffusion-related mixing times between 
                           diffusion-encoding block 1 and 2, block 2 and 3, 
                           in a DDE, TDE, et c setting. In DDE-based FEXI, 
                           we have only mde_tm1. 

- mde_b1, mde_b2, ...:     b-value per block in a DDE, TDE et c sequence
- mde_bt1, mde_bt2, ...:   b-tensor per block in a DDE, TDE et c sequence


Optional, may be calculated automatically:

- mde_tm12_ind:            Index according to mixing times
- mde_b1_ind, mde_b2_ind:  Index according to b1-values, b2-values, ...


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

