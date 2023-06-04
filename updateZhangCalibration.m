function [imgData, K] = updateZhangCalibration(imgData)
% INPUTS:
% imgData --> All the info about ALL the input imgs
% OUTPUTS:
% imgData --> All the info about ALL the input imgs with the updated H, P,
%             R and t for each img
% % K --> Normalized matrix of intrinsic parameters after radial
%         compensation
% OSS:
% 1) To use only inside Kevin_Marzio_project_*
% 2) Now we are using the refined_uv as the detected uv

    n_of_imgs = length(imgData);
   
    % Updating Homographies
    for ii=1:n_of_imgs  % 1 Homography matrx for each img
    
        XYpixel = imgData(ii).XYpixel;
    
        hom_col = ones(length(imgData(ii).XYmm), 1);
        hom_XYmm = [imgData(ii).XYmm, hom_col];
    
        imgData(ii).H = estimateH(XYpixel, hom_XYmm);
    end

    % Updating intrinsich parameter
    K = estimateIntrinsicStd(imgData);

    % Updating extrinsic parametrs e computing P
    for ii=1:n_of_imgs
        H = imgData(ii).H;
    
        [P, R, t] = estimateExtrinsicFromHK(H, K);
    
        imgData(ii).t = t;
        imgData(ii).R = R;
        imgData(ii).P = P;
    end
end