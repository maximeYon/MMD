Constrained diffusional variance decomposition (CODIVIDE)

The model features three components with Gaussian diffusion. The first
component is referred to as the free water component (FW), is isotropic
and has a diffusivity of 3 um2/ms. The second and third represent tissue
water, and has none and complete anisotropy, respectively, thereby having
values of d_delta fixed to 0 and 1. They share the same mean diffusivity,
which is fitted as a free parameter. 

In all, the model has four parameters: S0, f_FW, f_AT, D_I;T, representing
the non-diffusion weighted signal, the free-water signal fraction, the 
tissue-fraction of water with anisotropic diffusion, and the isotropic
diffusivity of the tissue component.

The model first appeared in Lampinen et al (2016) NeuroImage (minor revision)

