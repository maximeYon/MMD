function rout = fVectorNormalize(rin)

norm = sqrt(rin.x.^2 + rin.y.^2 + rin.z.^2);
eq.x = 'rout.x = rin.x./norm;'; Evaleq_xyz
