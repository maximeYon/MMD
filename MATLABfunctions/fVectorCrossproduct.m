function rout = fVectorCrossproduct(rin1,rin2)

rout.x = rin1.y.*rin2.z - rin1.z.*rin2.y;
rout.y = rin1.z.*rin2.x - rin1.x.*rin2.z;
rout.z = rin1.x.*rin2.y - rin1.y.*rin2.x;
