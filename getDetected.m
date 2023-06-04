function detected = getDetected(imgData_elem)
% INPUTS:
% imgData_elem --> All the info about a single img
% OUTPUTS:
% detected --> The pixels detected with detectCheckerboardPoints in the
%              input img (size is [N,2])
% OSS:
% 1) To use only inside Kevin_Marzio_project_*

    detected = imgData_elem.detectedXYpixel;
end