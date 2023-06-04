%% Intro

clc;
clear all;
close all;

n_of_imgs = 20;
% size_of_imgs = [480, 640];
checkerboard_size = [13, 14];
squaresize = 30;  % [mm]
visualize = false;

%% Loading imgs and detecting keypoints

% detect keypoints in the img
for ii=1:n_of_imgs
    imageFileName = fullfile('images',['image' num2str(ii) '.tif']);
  
    % we are defining an array of structures
    imgData(ii).name = imageFileName;
    imgData(ii).I = imread(imageFileName);
    imgData(ii).XYpixel = detectCheckerboardPoints(imgData(ii).I,...
                                                   'HighDistortion', false,...
                                                   'PartialDetections', true);
    badResult = sum(isnan(imgData(ii).XYpixel), 'all');
    if badResult
        disp([num2str(ii) ' HAS A HUGE PROBL']);
        imgData(ii).badResult = true;
        
    else
        disp([num2str(ii) ' is ok']);
        imgData(ii).badResult = false;
    end

    % to visualize
    if visualize
        hnd = figure('Name', imageFileName, 'NumberTitle', 'off');
        imshow(imgData(ii).I);
        hold on
      
        for jj=1:size(imgData(ii).XYpixel,1)
            % detected keypoints in the actual img
            x = imgData(ii).XYpixel(jj,1);
            y = imgData(ii).XYpixel(jj,2);
            
            plot(x, y, 'or')
            hndtxt=text(x, y, num2str(jj));
            set(hndtxt,'fontsize',12,'color','green');
        end
    
        pause
        close(hnd)
    end
end

% Keep only well detected imgs
idx = [imgData.badResult];
imgData = imgData(~idx);
n_of_imgs = sum(~idx);

% Storing cohordinates of each keypoint in the checkerboard plane
for ii = 1:n_of_imgs
    
    XYpixel = imgData(ii).XYpixel;  % img keypoints 
    clear Xmm Ymm;
    
    if visualize
        hnd = figure('Name', imgData(ii).name, 'NumberTitle', 'off');
        imshow(imgData(ii).I)
    end
    
    for jj=1:size(XYpixel,1)  % length(XYpixel)
        
        [row, col] = ind2sub([12, 13], jj);  % linear indexes
        Xmm = (col-1)*squaresize;
        Ymm = (row-1)*squaresize;
        imgData(ii).XYmm(jj,:) = [Xmm, Ymm];  % riga jj, tutta la colonna
       
        if visualize
            hndtext = text(XYpixel(jj,1), XYpixel(jj,2), num2str([Xmm, Ymm]));
            set(hndtext, 'fontsize', 12, 'color', 'cyan');
        end
    end

    if visualize
        pause
        close(hnd)
    end

end

% save('keypoints_function');

%% Estimating homography

% load('keypoints_function');

for ii=1:n_of_imgs  % 1 Homography matrx for each img

    XYpixel = imgData(ii).XYpixel;

    hom_col = ones(length(imgData(ii).XYmm), 1);
    hom_XYmm = [imgData(ii).XYmm, hom_col];

    imgData(ii).H = estimateH(XYpixel, hom_XYmm);
end

%% Finding intrinsic parameters (CHOLESKY STANDARD)

K = estimateIntrinsicStd(imgData);

%% Finding extrinsic parametrs e computing P

for ii=1:n_of_imgs
    H = imgData(ii).H;

    [P, R, t] = estimateExtrinsicFromHK(H, K);

    imgData(ii).t = t;
    imgData(ii).R = R;
    imgData(ii).P = P;

end

% save('keypoints_function_2')

%% Reprojection error...

% load('keypoints_function_2')

% First estimate is assumed ideal without radial distorsion
k1 = 0;     k2 = 0;

for ii=1:n_of_imgs
    % Secures detected points
    imgData(ii).detectedXYpixel = imgData(ii).XYpixel;
    
    % Secures original reprojection error (in absence of radial dist. compensation)
    imgData(ii).eps = computeReprojError(imgData(ii),...  % on which img
                                         K,...            % Intrinsic param matrix
                                         k1, k2);         % Radial distorsion params
end

%% ... and detail of an img

kk = randi(n_of_imgs);
kk = 15;

disp(['Reprojection error, without radial dist compensation for img', ...
      imgData(kk).name, ' is:', ...
      num2str(imgData(kk).eps)]);

hnd = figure('Name', imgData(kk).name, 'NumberTitle', 'off');
imshow(imgData(kk).I);
% impixelinfo
hold on

detected_points = getDetected(imgData(kk));
projected_points = getProjected(imgData(kk));

% figure with detected and projected points
for jj=1:length(detected_points)

    % measured pixel cohordinates
    u = detected_points(jj, 1);
    v = detected_points(jj, 2);
    
    % estimated pixel cohordinates (first estimate is assumed ideal)
    est_u = projected_points(jj,1);
    est_v = projected_points(jj,2);

    % visualization
    plot(u, v, 'xg')
    plot(est_u, est_v, 'xb')
    legend('detected', 'estimated from P')

end

% pause;
% close(hnd);
hold off

%% Radial distorsion compensation

% load('keypoints_function_2')

compensation_step = 1;
mean_eps_array = [];

% -------------------------------------------------------------------------
% 2) first estimate of k1, k2 (radial distorsion parameters) (L2, pag 69)
% -------------------------------------------------------------------------

% distorted cohordinates (obtained from the img)
all_detected_points = getAllDetected(imgData);

% undistorted cohordinates (assuming P "perfect")
all_projected_points = getAllProjected(imgData);

[k1, k2] = estimateRadDistParams(K, all_projected_points, all_detected_points);


% -------------------------------------------------------------------------
% 3) compensation of radial distorsion                       (L2, pag 69)
% -------------------------------------------------------------------------

% Refine projected points
for ii=1:n_of_imgs
    imgData(ii).XYpixel = compensateRadialDistorsion(imgData(ii), K, k1, k2);
end

% -------------------------------------------------------------------------
% Iterate radial compensation until convergence...
% -------------------------------------------------------------------------

trheshold = 1;

mean_eps = meanReprErr(imgData, K, k1, k2);
mean_eps_array = [mean_eps_array, mean_eps];

compensation_step = 2;


while true

    mean_eps_old = mean_eps;

    % Update the result obtained from Zhang's calibration
    [imgData, K] = updateZhangCalibration(imgData);
    
    % computing the new radial distorsion parameters...
    all_detected_points = getAllDetected(imgData);
    all_projected_points = getAllProjected(imgData);
    [k1, k2] = estimateRadDistParams(K, all_projected_points, all_detected_points);
    
    % refining projected points...
    for ii=1:n_of_imgs
        imgData(ii).XYpixel = compensateRadialDistorsion(imgData(ii), K, k1, k2);
    end
    
    % checking improvements
    mean_eps = meanReprErr(imgData, K, k1, k2);
    mean_eps_array = [mean_eps_array, mean_eps];
    d_eps = mean_eps-mean_eps_old;

    if(abs(d_eps) < trheshold)
        disp(['Stopped ad step ', num2str(compensation_step)])
        disp(['Results are k1=', num2str(k1), ', k2=', num2str(k2)]);
        disp(['mean_eps is ', num2str(mean_eps)]);
        break
    else
        disp(['At step ', num2str(compensation_step),...
              ' results are: d_eps=', num2str(d_eps), ', mean_eps=', num2str(mean_eps)]);
        disp(['At step ', num2str(compensation_step),...
              ' in ', imgData(kk).name,  ' sqrt(eps)=', num2str(sqrt(computeReprojError(imgData(kk), K, k1, k2)))]);
    end

    compensation_step = compensation_step + 1;

end

% save('keypoints_function_3')

%% Total reprojection error comparison

for ii=1:n_of_imgs
    % Reprojection error (after radial dist. compensation)
    imgData(ii).eps_new = computeReprojError(imgData(ii),...
                                             K,...
                                             k1, k2);

    % Notify user if something bad occured
    if (imgData(ii).eps_new > imgData(ii).eps)
        disp(['For ', imgData(ii).name, ' reproj. error is INCREASED!!!']);
    end
end

total_eps = 0;
total_eps_new = 0;

% total reproj. error
for ii=1:n_of_imgs
    total_eps = total_eps + imgData(ii).eps;
    total_eps_new = total_eps_new + imgData(ii).eps_new;
end

% mean reproj. error
mean_eps = total_eps/n_of_imgs;    mean_eps_new = total_eps_new/n_of_imgs;

disp(['Mean reprojection error before radial compensation is: ',...
      num2str(mean_eps)]);
disp(['For Image 17 was: ', num2str(imgData(15).eps)]);
disp(['Mean reprojection error after radial compensation is: ',...
      num2str(mean_eps_new)]);
disp(['The difference is: ', num2str(mean_eps_new - mean_eps)]);
disp(['For Image 17 is: ', num2str(imgData(15).eps_new)]);

% plot of reprojection error during compensation procedure
hnd = figure('Name', 'Mean reprojection error', 'NumberTitle', 'off');
n = 1:compensation_step;
plot(n, mean_eps_array)
xlabel('Compensation step')
ylabel('Mean reprojection error')


%% figure of detected and projected points after compensation

hnd = figure('Name', [imgData(kk).name, ' after readial compensation'], 'NumberTitle', 'off');
imshow(imgData(kk).I);
% impixelinfo
hold on

alpha_u = K(1,1);
alpha_v = K(2,2);
u0 = K(1,3);
v0 = K(2,3);

detected_points = imgData(kk).detectedXYpixel;
projected_points = getProjected(imgData(kk));

% figure with detected and projected points
for jj=1:length(detected_points)

    % measured pixel cohordinates
    u = detected_points(jj, 1);
    v = detected_points(jj, 2);
    
    % estimated pixel cohordinates (P)
    est_u_id = projected_points(jj,1);
    est_v_id = projected_points(jj,2);
    
    % estimated pixel cohordinates (P + radial distorsion)
    rd2 = ( (est_u_id-u0)/alpha_u )^2 + ( (est_v_id-v0)/alpha_v )^2;
    rd4 = rd2^2;
    est_u_hat = (est_u_id - u0)*(1 + k1*rd2 + k2*rd4) + u0;
    est_v_hat = (est_v_id - v0)*(1 + k1*rd2 + k2*rd4) + v0;

    % visualization
    plot(u, v, 'xg')
    plot(est_u_id, est_v_id, 'xb')
    plot(est_u_hat, est_v_hat, 'xr')

    legend('detected', 'P', 'P + rad. dist.')
end

% pause;
% close(hnd);
hold off

%% Superimpose an object on imgs

for kk=1:n_of_imgs
    P = imgData(kk).P;
    
    hnd = figure('Name', [imgData(kk).name, ' with superimposition'], 'NumberTitle', 'off');
    imshow(imgData(kk).I);
    hold on
    
    % where to place the cylinder (center of the checkerboard)
    C = (checkerboard_size - 1).*squaresize./2;
    
    % cylinder properties
    vtheta = linspace(0,2*pi);  %0:0.01:2*pi
    C_x = C(1);    C_y = C(2);    r = 75;
    vx = C_x + r*cos(vtheta);      % Xmm
    vy = C_y + r*sin(vtheta);      % Ymm
    
    % base of the cylinder (in the checkerboard plane)
    vz = zeros(1,length(vtheta));  % Zmm
    hom_row = ones(1,length(vtheta));
    hom_XYZmm = [vx; vy; -vz; hom_row];
    hom_XYpixel = P * hom_XYZmm;
    XYpixel = [hom_XYpixel(1,:)./hom_XYpixel(3,:);...
               hom_XYpixel(2,:)./hom_XYpixel(3,:)];
    
    % ... and its colour specs.
    f1=fill(XYpixel(1,:),XYpixel(2,:),'r');
    set(f1,'facealpha',.5)
    
    % roof of the cylinder...
    vz = 3*squaresize*ones(1,length(vtheta));  % Zmm
    hom_XYZmm = [vx; vy; -vz; hom_row];
    hom_XYpixel = P * hom_XYZmm;
    XYpixel = [hom_XYpixel(1,:)./hom_XYpixel(3,:);...
               hom_XYpixel(2,:)./hom_XYpixel(3,:)];
    
    % ... and its colour specs.
    f2=fill(XYpixel(1,:),XYpixel(2,:),'g');
    set(f2,'facealpha',.5)

    % Creating the axes...
    v_ax = 0:0.01:2*squaresize;
    x_ax = 0 + v_ax;
    y_ax = 0 + v_ax;
    z_ax = 0 + v_ax;
    hom_row = ones(1,length(v_ax));
    O_row = zeros(1, length(v_ax));
    
    % in real-world hom-coordinates...
    hom_x_ax_mm = [v_ax; O_row; O_row; hom_row];
    hom_y_ax_mm = [O_row; v_ax; O_row; hom_row];
    hom_z_ax_mm = [O_row; O_row; v_ax; hom_row];

    % in img hom-coordinates...
    hom_x_ax_pixel = P * hom_x_ax_mm;
    hom_y_ax_pixel = P * hom_y_ax_mm;
    hom_z_ax_pixel = P * hom_z_ax_mm;

    % in img coordinates...
    x_ax_pixel = [hom_x_ax_pixel(1,:)./hom_x_ax_pixel(3,:);...
               hom_x_ax_pixel(2,:)./hom_x_ax_pixel(3,:)];
    y_ax_pixel = [hom_y_ax_pixel(1,:)./hom_y_ax_pixel(3,:);...
               hom_y_ax_pixel(2,:)./hom_y_ax_pixel(3,:)];
    z_ax_pixel = [hom_z_ax_pixel(1,:)./hom_z_ax_pixel(3,:);...
               hom_z_ax_pixel(2,:)./hom_z_ax_pixel(3,:)];

    % and plot them
    plot(x_ax_pixel(1,:), x_ax_pixel(2,:), 'r', LineWidth=1);
    text(x_ax_pixel(1,end), x_ax_pixel(2,end), 'X');
    plot(y_ax_pixel(1,:), y_ax_pixel(2,:), 'g', LineWidth=1);
    text(y_ax_pixel(1,end), y_ax_pixel(2,end), 'Y');
    plot(z_ax_pixel(1,:), z_ax_pixel(2,:), 'b', LineWidth=1);
    text(z_ax_pixel(1,end), z_ax_pixel(2,end), 'Z');

    hold off
end