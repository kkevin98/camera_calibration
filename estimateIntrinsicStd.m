function K = estimateIntrinsicStd(imgData)
% OUTPUTS:
% K --> Estimated matrix of intrinsic parameters of the camera
% OSS:
% 1) To use only inside Kevin_Marzio_project_*
    
    % "Creating" the system (Vb=o), to find B...
    n_of_imgs = length(imgData);
    V = zeros(2*n_of_imgs, 6);

    for ii=1:n_of_imgs
        H = imgData(ii).H;
    
        % creating v_ij
        v11 = [     H(1,1)*H(1,1);
               H(1,1)*H(2,1) + H(2,1)*H(1,1);
                    H(2,1)*H(2,1);
               H(3,1)*H(1,1) + H(1,1)*H(3,1);
               H(3,1)*H(2,1) + H(2,1)*H(3,1);
                    H(3,1)*H(3,1)];
        v12 = [     H(1,1)*H(1,2);
               H(1,1)*H(2,2) + H(2,1)*H(1,2);
                    H(2,1)*H(2,2);
               H(3,1)*H(1,2) + H(1,1)*H(3,2);
               H(3,1)*H(2,2) + H(2,1)*H(3,2);
                    H(3,1)*H(3,2)];
        v22 = [     H(1,2)*H(1,2);
               H(1,2)*H(2,2) + H(2,2)*H(1,2);
                    H(2,2)*H(2,2);
               H(3,2)*H(1,2) + H(1,2)*H(3,2);
               H(3,2)*H(2,2) + H(2,2)*H(3,2);
                    H(3,2)*H(3,2)];
    
        % ... stacking eqs. for each plane/img...
        V(2*ii-1,:) = v12';
        V(2*ii,:) = (v11 - v22)';
    end
    
    % ... properly solve the system to find B
    [~, ~, S] = svd(V);
    b = S(:, end);
    B = [b(1), b(2), b(4);...
         b(2), b(3), b(5);...
         b(4), b(5), b(6)];
    
    % Find K
    [L, not_pos_def] = chol(B, "lower");
    if not_pos_def
        %disp(['B is not positive definite']);
        L = chol(-B, "lower");
    else
        %disp(['B is positive definite']);
    end
    K = inv(L');
    K = 1/K(3,3).*K;  % normalize
end