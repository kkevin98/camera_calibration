function H = estimateH(XYpixel, hom_XYmm)
% INPUTS:
% XYpixel --> size(XYpixel_norm) should be [N,2]
% XYmm --> size(hom_XYmm) should be [N,3]
% OUTPUTS:
% H --> Estimated Homography
    
    % Creating the system (A*h=0) to find H...
    n_of_corr = size(XYpixel,1);
    A = zeros(2*n_of_corr, 9);
    O = [0; 0; 0];

    for ii=1:n_of_corr
        u = XYpixel(ii, 1);
        v = XYpixel(ii, 2);
        m = hom_XYmm(ii,:)';
             
        %  ... stacking eqs. for each correspondance
        A(2*ii-1,:) = [m', O', -u*m'];
        A(2*ii,:) = [O', m', -v*m'];
    end
    
    % Properly solve the system (+ enforcing 8DOF in H)
    % [~, ~, V] = svd((A')*A);  % Collins proposal
    [~, ~, V] = svd(A);
    h = V(:, end);
    H = reshape(h, [3, 3])';

    % Making the position of the camera above and not below the plane
    if det(H)<0
        H = -H;
    end
end