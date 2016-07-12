function [L, v_1] = tm_1x6_eigvals(v_1x6, c_method)
% function [L, v_1] = tm_1x6_eigvals(v_1x6, c_method)
%
%
% c_method
%   1 - Cardano's fast method
%   2 - eigs in Matlab
if (nargin < 2), c_method = 1; end

switch (c_method)
    
    case 1
        % Use Cardano's method to diagonalization
        c_0 = ...
            v_1x6(:,1) .* abs(v_1x6(:,6)).^2/2 + ...
            v_1x6(:,2) .* abs(v_1x6(:,5)).^2/2 + ...
            v_1x6(:,3) .* abs(v_1x6(:,4)).^2/2 - ...
            v_1x6(:,1) .* v_1x6(:,2) .* v_1x6(:,3) - ...
            2 * real(conj(v_1x6(:,5)) .* v_1x6(:,6) .* v_1x6(:,4)) / sqrt(2)^3;
        
        c_1 = v_1x6(:,1).*v_1x6(:,2) + v_1x6(:,1).*v_1x6(:,3) + v_1x6(:,2).*v_1x6(:,3) - ...
            sum(abs(v_1x6(:,4:6)).^2/2,2);
        
        c_2 = -sum(v_1x6(:,1:3),2);
        
        
        p = c_2.^2 - 3 * c_1;
        q = -27/2 * c_0 - c_2.^3  + 9/2 * c_2 .* c_1;
        
        phi = 1/3 * atan( sqrt( 27 * (1/4 * c_1.^2 .* (p - c_1) + c_0 .* (q + 27/4 * c_0))) ./ q);
        
        phi = phi + pi * (q < 0);
        
        
        x1 = 2 * cos(phi);
        x2 = -cos(phi) - sqrt(3) * sin(phi);
        x3 = -cos(phi) + sqrt(3) * sin(phi);
        
        ps = sqrt(p) / 3;
        
        L = [ps .* x1 - 1/3 * c_2 ps .* x2 - 1/3 * c_2 ps .* x3 - 1/3 * c_2] ;
        
        L = sort(L, 2, 'descend');
        
        % if desired, also compute eigenvectors
        %   fast code from CF
        if (nargout > 1)
            
            %             t11 = v_1x6(:,1);
            %             t22 = v_1x6(:,2);
            %             t33 = v_1x6(:,3);
            %
            %             t12 = v_1x6(:,4) / sqrt(2);
            %             t13 = v_1x6(:,5) / sqrt(2);
            %             t23 = v_1x6(:,6) / sqrt(2);
            %
            %             l2 = L(:,2);
            %             l3 = L(:,3);
            %
            %             x =  (t11-l2).*(t11-l3)+t12.^2+t13.^2;
            %             y =  t12.*(t11-l3)+(t22-l2).*t12+t13.*t23;
            %             z =  t13.*(t11-l3)+t12.*t23+(t33-l2).*t13;
            
            x =  (v_1x6(:,1)-L(:,2)).*(v_1x6(:,1)-L(:,3))+v_1x6(:,4).^2/2+v_1x6(:,5).^2/2;
            y =  v_1x6(:,4).*(v_1x6(:,1)-L(:,3))/sqrt(2)+(v_1x6(:,2)-L(:,2)).*v_1x6(:,4)/sqrt(2)+v_1x6(:,5).*v_1x6(:,6)/2;
            z =  v_1x6(:,5)/sqrt(2).*(v_1x6(:,1)-L(:,3))+v_1x6(:,4)/sqrt(2).*v_1x6(:,6)/sqrt(2)+(v_1x6(:,3)-L(:,2)).*v_1x6(:,5)/sqrt(2);
            
            tnorm = 1./sqrt(x.*x + y.*y + z.*z +eps);
            e1x = tnorm.*x;
            e1y = tnorm.*y;
            e1z = tnorm.*z;
            
            v_1 = [e1x e1y e1z];
            
            
            % x =  (t11-l1).*(t11-l3)+t12.^2+t13.^2;
            % y =  t12.*(t11-l3)+(t22-l1).*t12+t13.*t23;
            % z =  t13.*(t11-l3)+t12.*t23+(t33-l1).*t13;
            % tnorm = 1./sqrt(x.*x + y.*y + z.*z +eps);
            % e2(:,:,:,1) = tnorm.*x;
            % e2(:,:,:,2) = tnorm.*y;
            % e2(:,:,:,3) = tnorm.*z;
            %
            % x =  (t11-l1).*(t11-l2)+t12.^2+t13.^2;
            % y =  t12.*(t11-l2)+(t22-l1).*t12+t13.*t23;
            % z =  t13.*(t11-l2)+t12.*t23+(t33-l1).*t13;
            % tnorm = 1./sqrt(x.*x + y.*y + z.*z +eps);
            % e3(:,:,:,1) = tnorm.*x;
            % e3(:,:,:,2) = tnorm.*y;
            % e3(:,:,:,3) = tnorm.*z;
            
        end
        
    case 2
        
        L = zeros(size(v_1x6,1), 3);
        v_1    = zeros(size(v_1x6,1), 3);
        
        for c = 1:size(v_1x6,1)
            t_3x3 = tm_1x6_to_3x3(v_1x6(c,:));
            [a,b] = eigs(t_3x3);
            
            [~,ind] = sort(diag(b), 'descend');
            a = a(:,ind);
            b = b(ind,ind);
            
            L(c,:) = diag(b);
            
            [~,ind] = max(diag(b));
            v_1(c,:) = a(:,ind);
            
        end
end
