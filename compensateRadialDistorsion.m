function refined_projection = compensateRadialDistorsion(imgData_elem, K, k1, k2)
% INPUTS:
% imgData_elem --> All the info about a single img
% K --> Normalized matrix of instrinsic parameters
% k1 --> 1st parameter of radial distorsion model
% k2 --> 2nd parameter of radial distorsion model
% OUTPUTS:
% refined_projection --> Refined version of the pixels projected in the 
%                        input img starting from the real world points. 
%                        This points does not take into account the radial 
%                        distorsion effects (size is [N,2])
% OSS:
% 1) To use only inside Kevin_Marzio_project_*
    
    detected_uv = getDetected(imgData_elem);        % uv_hat
    projected_points = getProjected(imgData_elem);  % old uv
    imgDim = size(imgData_elem.I);                  % [height, width]
    n_of_points = length(detected_uv);

    alpha_u = K(1,1);
    alpha_v = K(2,2);
    u0 = K(1,3);
    v0 = K(2,3);
    
    refined_projection = zeros(n_of_points, 2);

    for jj=1:n_of_points
        
        % distorted cohordinates
        x_hat = (detected_uv(jj,1) - u0)/alpha_u;
        y_hat = (detected_uv(jj,2) - v0)/alpha_v;

        % inital guess to solve the system
        x0 = (projected_points(jj,1) - u0)/alpha_u;
        y0 = (projected_points(jj,2) - v0)/alpha_v;

        % sytem to be solved to find the "ideal" xy (L2 pag 69)
        non_linear_syst = @(xy)radial_dist_compens_syst(xy,...
                                                        x_hat, y_hat,...
                                                        k1, k2);
        
        % the solution
        options = optimset('Display','off');
        xy = fsolve(non_linear_syst, [x0; y0], options);  % normalized cohordinates
        u = alpha_u*xy(1)+u0;                             % standard cohordinates
        v = alpha_v*xy(2)+v0;                             %     standard cohordinates

        % checking the results and notify the user if some cohordinates
        % does not belong to the img
        if (u>0) && (u<imgDim(2)) && (v>0) && (v<imgDim(1))
            refined_projection(jj,:) = [u, v];
        else
            disp(['A problem in ' imgData_elem.name ' occured:'...
                  ' refined cohordinates does not lie in the img!!!'])
            disp(['Obtained [u,v] are: [', num2str(u), ', ' num2str(v), ']'])
            refined_projection(jj,:) = [u, v];
        end
    end

    
    function S = radial_dist_compens_syst(xy, x_hat, y_hat, k1, k2)
        x = xy(1);
        y = xy(2);
    
        S(1) = x*(1 + k1*(x^2+y^2) + k2*(x^4+2*(x^2)*(y^2)+y^4) ) - x_hat;
        S(2) = y*(1 + k1*(x^2+y^2) + k2*(x^4+2*(x^2)*(y^2)+y^4) ) - y_hat;
    end
end