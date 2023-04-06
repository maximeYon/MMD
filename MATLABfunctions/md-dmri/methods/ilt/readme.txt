Neurite density imaging (NDI) model

This model is a simplification of the neurite density and dispersion 
imaging model by Zhang et al (2012) NeuroImage. The simplification 
comes by utilizing the powder-averaged signal, which removes the dependency
of the signal on orientation and dispersion. Thus, only three parameters
remain to be fitted: S0, f_CSF, and f_IC, corresponding to the non-diffusion
weighted signal, the free-water fraction, and the neurite density. Note that
on tissues that do not fulfil the tortuosity assumptoin, such as in gray
matter or in neurites, the neurite fraction will be biased. Read more
about this in Lampinen et al (2016) NeuroImage (submitted/revised). 

Note that this implementation is slightly different from Lampinen's. 

