function coord_out = fPutMoleculesInSameFrameCoords(coord_in,frame_in,coord_ref,frame_ref,box)

coord_out = coord_in;

tol = 1e-3;

Natoms = numel(coord_ref.x(1,:));
for natom = 1:Natoms

    if abs((coord_out.x(frame_in,natom) - coord_ref.x(frame_ref,natom))) > (1-tol)*box.xx(1,frame_in)
        multiplier = (coord_out.x(frame_in,natom) - coord_ref.x(frame_ref,natom))/box.xx(1,frame_in);
        coord_out.x(:,natom) = coord_in.x(:,natom) - multiplier*box.xx(1,frame_in);
    end
    if abs((coord_out.y(frame_in,natom) - coord_ref.y(frame_ref,natom))) > (1-tol)*box.yy(1,frame_in)
        multiplier = (coord_out.y(frame_in,natom) - coord_ref.y(frame_ref,natom))/box.yy(1,frame_in);
        coord_out.y(:,natom) = coord_in.y(:,natom) - multiplier*box.yy(1,frame_in);
    end
    if abs((coord_out.z(frame_in,natom) - coord_ref.z(frame_ref,natom))) > (1-tol)*box.zz(1,frame_in)
        multiplier = (coord_out.z(frame_in,natom) - coord_ref.z(frame_ref,natom))/box.zz(1,frame_in);
        coord_out.z(:,natom) = coord_in.z(:,natom) - multiplier*box.zz(1,frame_in);
    end

end
