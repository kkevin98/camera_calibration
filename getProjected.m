function projected = getProjected(imgData_elem)
% INPUTS:
% imgData_elem --> All the info about a single img
% OUTPUTS:
% projected --> The pixels projected in the input img starting from the 
%               real world points. This points does not take into account
%               the radial distorsion effects (size is [N,2])
% OSS:
% 1) To use only inside Kevin_Marzio_project_*
    
    XYmm = imgData_elem.XYmm;
    n_of_points = length(XYmm);
    P = imgData_elem.P;
    p1t = P(1,:);    p2t = P(2,:);    p3t = P(3,:);  % row of P
    projected = zeros(n_of_points, 2);

    % computing projection
    for jj=1:n_of_points
    
        % homogeneous world points
        m_x = XYmm(jj, 1);
        m_y = XYmm(jj, 2);
        m_z = 0;
        m = [m_x; m_y; m_z; 1];
        
        % estimated pixel cohordinates
        projected(jj, 1) = (p1t*m) / (p3t*m);
        projected(jj, 2) = (p2t*m) / (p3t*m);
   
    end
end