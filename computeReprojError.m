function eps = computeReprojError(imgData_elem, K, k1, k2)
% INPUTS:
% imgData_elem --> All the info about a single img
% K --> Normalized matrix of intrinsic parameters
% k1 --> 1st parameter of radial distorsion model
% k2 --> 2nd parameter of radial distorsion model
% OUTPUTS:
% eps --> reprojection error (L2, pag 62)
% OSS:
% 1) To use only inside Kevin_Marzio_project_*

    n_of_points = length(imgData_elem.XYmm);
    alpha_u = K(1,1);
    alpha_v = K(2,2);
    u0 = K(1,3);
    v0 = K(2,3);

    % If radial distrsion is not present P would project the points in...
    ideal_estimated_img_points = getProjected(imgData_elem);
    
    % Real projections in the img due to radial ditstorsion...
    estimated_img_points = zeros(n_of_points, 2);  % preallocating

    for ii=1:n_of_points
        u_id = ideal_estimated_img_points(ii,1);
        v_id = ideal_estimated_img_points(ii,2);

        rd2 = ( (u_id-u0)/alpha_u )^2 + ( (v_id-v0)/alpha_v )^2;
        rd4 = rd2^2;

        u_hat = (u_id - u0)*(1 + k1*rd2 + k2*rd4) + u0;
        v_hat = (v_id - v0)*(1 + k1*rd2 + k2*rd4) + v0;

        estimated_img_points(ii,1) = u_hat;
        estimated_img_points(ii,2) = v_hat;
    end

    % Detected points in the img are...
    detected_img_points = getDetected(imgData_elem);

    % Reprojection error as in L2, pag 62
    squared_errors = (estimated_img_points - detected_img_points).^2;
    eps = sum(squared_errors, "all");
end