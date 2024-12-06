function gro_out = fPutMoleculesInBoxGro(gro_in)

gro = gro_in;

molnr = unique(gro.molnr);
for ncount = 1:numel(molnr)
    index = find(gro.molnr==molnr(ncount));

    thresh = .5;
    for m = 2:length(index)
        if (gro.x(index(m))-gro.x(index(m-1))) > thresh*gro.boxxx
            gro.x(index(m)) = gro.x(index(m)) - gro.boxxx;
        end
        if (gro.x(index(m))-gro.x(index(m-1))) < -thresh*gro.boxxx
            gro.x(index(m)) = gro.x(index(m)) + gro.boxxx;
        end
        if (gro.y(index(m))-gro.y(index(m-1))) > thresh*gro.boxyy
            gro.y(index(m)) = gro.y(index(m)) - gro.boxyy;
        end
        if (gro.y(index(m))-gro.y(index(m-1))) < -thresh*gro.boxyy
            gro.y(index(m)) = gro.y(index(m)) + gro.boxyy;
        end
        if (gro.z(index(m))-gro.z(index(m-1))) > thresh*gro.boxzz
            gro.z(index(m)) = gro.z(index(m)) - gro.boxzz;
        end
        if (gro.z(index(m))-gro.z(index(m-1))) < -thresh*gro.boxzz
            gro.z(index(m)) = gro.z(index(m)) + gro.boxzz;
        end
    end

    if mean(gro.x(index)) > gro.boxxx
        gro.x(index) = gro.x(index) - gro.boxxx;
    end
    if mean(gro.x(index)) < 0
        gro.x(index) = gro.x(index) + gro.boxxx;
    end
    if mean(gro.y(index)) > gro.boxyy
        gro.y(index) = gro.y(index) - gro.boxyy;
    end
    if mean(gro.y(index)) < 0
        gro.y(index) = gro.y(index) + gro.boxyy;
    end
    if mean(gro.z(index)) > gro.boxzz
        gro.z(index) = gro.z(index) - gro.boxzz;
    end
    if mean(gro.z(index)) < 0
        gro.z(index) = gro.z(index) + gro.boxzz;
    end

end

gro_out = gro;