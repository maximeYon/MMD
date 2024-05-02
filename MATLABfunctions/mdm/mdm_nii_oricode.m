function o = mdm_nii_oricode(h)
% function o = mdm_nii_oricode(h)
%
% Find the principal orientation of the nifti header

if (h.sform_code > 0)
    t = [h.srow_x(1:3) h.srow_y(1:3) h.srow_z(1:3)]';
elseif (h.qform_code > 0)
    b = h.quatern_b;
    c = h.quatern_c;
    d = h.quatern_d;
    
    q = sum([b c d].^2);
    if (1 - q < 0)
        if (abs(1-q) < 1e-5)
            a = 0;
        else
            error('Incorrect quaternions');
        end
    else
        a = sqrt(1-q);
    end
    
    qfac = h.pixdim(1);    
    if (qfac == 0), qfac = 1; end
    
    i = h.pixdim(2);
    j = h.pixdim(3);
    k = qfac * h.pixdim(4);
    
    R = [...
        a*a+b*b-c*c-d*d     2*b*c-2*a*d        2*b*d+2*a*c
        2*b*c+2*a*d         a*a+c*c-b*b-d*d    2*c*d-2*a*b
        2*b*d-2*a*c         2*c*d+2*a*b        a*a+d*d-c*c-b*b];
    
    if (det(R) == 0)
        tolerance = 1;
        R_sort = sort(abs(R(:)));
        R( ( abs(R) < tolerance*min(R_sort(end-2:end)) ) ) = 0;
    end
    
    t = R * diag([i j k]);
    
else % assume it is in LAS form
    t = diag([-1 1 1]); %'LAS';
end

p = {'RAS', 'LPI'}; % right handed systems

o = [];
for i = 1:3
    [~,i_o] = max(abs(t(:,i)));
    i_s = (t(i_o, i) < 0) + 1;
    o = [o p{i_s}(i_o)];
end
