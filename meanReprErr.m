function eps_mean = meanReprErr(imgData, K, k1, k2)
% INPUTS:
% imgData_elem --> All the info about ALL the input imgs
% K --> Normalized matrix of intrinsic parameters
% k1 --> 1st parameter of radial distorsion model
% k2 --> 2nd parameter of radial distorsion model
% OUTPUTS:
% eps_mean --> mean reprojection error computed 
%              between all imgs(L2, pag 62)
% OSS:
% 1) To use only inside Kevin_Marzio_project_*

    n_of_imgs = length(imgData);
    eps_total = 0;

    for ii=1:n_of_imgs
        eps_total = eps_total + computeReprojError(imgData(ii), K, k1, k2);
    end
    eps_mean = eps_total/n_of_imgs;
end