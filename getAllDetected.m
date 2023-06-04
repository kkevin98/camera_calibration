function all_detected = getAllDetected(imgData)
% INPUTS:
% imgData --> All the info about ALL the input imgs
% OUTPUTS:
% all_detected --> The pixels detected with detectCheckerboardPoints in ALL
%                  the input imgs (size is [N,2])
% OSS:
% 1) To use only inside Kevin_Marzio_project_*
    
    n_of_imgs = length(imgData);
    all_detected = [];

    for ii=1:n_of_imgs
        all_detected = [all_detected;...
                        getDetected(imgData(ii))];
    end

end