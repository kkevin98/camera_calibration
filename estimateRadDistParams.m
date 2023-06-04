function [k1, k2]= estimateRadDistParams(K, uv_ideal, uv_hat)
% INPUTS:
% K --> Normalized matrix of instrinsic parameters
% uv_ideal --> undistorted pixel cohordinates (computed asusming P
%              "perfect").
%              size(uv_ideal) must be [N*2]
% uv_hat --> distorted pixel cohordinates (detected from the distorted img)
%            size(uv_hat) must be [N*2]
% OUTPUTS:
% k1 --> first param. related to radial distorsion
% k2 --> second param. related to radial distorsion

    % needed intrinsic params
    alpha_u = K(1,1);
    alpha_v = K(2,2);
    u0 = K(1,3);
    v0 = K(2,3);

    n_of_points = length(uv_ideal);
    A = zeros(2*n_of_points, 2);
    b = zeros(2*n_of_points, 1);
    
    % Creating the system A*[k1; k2]=b ...
    for ii=1:n_of_points
        % undistorted cohordinates
        u_ideal = uv_ideal(ii, 1);
        v_ideal = uv_ideal(ii, 2);

        % distorted cohordinates
        u_hat = uv_hat(ii, 1);
        v_hat = uv_hat(ii, 2);

        rd2 = ( (u_ideal - u0)/alpha_u )^2 +...
              ( (v_ideal - v0)/alpha_v )^2;
        
        % ... stacking the equations...
        A(2*ii-1,1)=(u_ideal-u0)*rd2;    A(2*ii-1,2)=(u_ideal-u0)*(rd2^2);
        A(2*ii,1) = (v_ideal-v0)*rd2;    A(2*ii,2)=(v_ideal-v0)*(rd2^2);
        
        b(2*ii-1, 1) = u_hat-u_ideal;
        b(2*ii, 1) = v_hat-v_ideal;
    end

    % ... solving the system
    rad_dist_params = ((A')*A) \ (A') * b;
    k1 = rad_dist_params(1);
    k2 = rad_dist_params(2);

end