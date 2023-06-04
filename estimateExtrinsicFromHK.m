function [P, R, t] = estimateExtrinsicFromHK(H, K)
% INPUTS:
% H --> Homography matrix
% K --> Normalized matrix of instrinsic parameters
% OUTPUTS:
% P --> Perspective projection matrix
% R --> Rotational matrix associated to P
% t --> traslational vector associated to P
    
    h1 = H(:,1);    h2 = H(:,2);    h3 = H(:,3);

    lambda1 = 1/norm(K\h1);
    lambda2 = 1/norm(K\h2);
    % rough approximation bcs lambda1 is not perfectly equal to lambda2
    lambda = (lambda1 + lambda2)/2;

    r1 = 1/lambda * K\h1;
    r2 = 1/lambda * K\h2;
    t = 1/lambda * K\h3;
    r3 = cross(r1, r2);
    
    % Orthogonal Procrustes problem to get orthogonal R
    R_est = [r1, r2, r3];
    [U, ~, V] = svd(R_est);
    R = U*V';
    
    % Warning user if we get a reflection instead of a rotation
    if det(R)>0
        % disp(['det(R) = ' num2str(det(R))]);
    else
        disp(['ATTENTION det(R) = ' num2str(det(R)) ' !!!']);
    end

    P = K*[R, t];
    
    % Warning user if  P is not a proj. matrix (Faugeras, L2 pag 43)
    Q = P([1:3],[1:3]);
    tol = 1.e-6; 
    if abs(det(Q)) < tol
        disp(['P is a non valid proj. matrix!!'])
    end
end