function all_projected = getAllProjected(imgData)
% INPUTS:
% imgData --> All the info about ALL the input imgs
% OUTPUTS:
% all_projected --> ALL the pixels projected in ALL the input imgs starting 
%                   from the real world points. This points does not take 
%                   into account the radial distorsion effects (size is [N,2])
% OSS:
% 1) To use only inside Kevin_Marzio_project_*
    
    n_of_imgs = length(imgData);
    all_projected = [];

    for ii=1:n_of_imgs

        % projected points contained on the current img
        this_projected = getProjected(imgData(ii));

        all_projected = [all_projected;...
                         this_projected];
    end

end