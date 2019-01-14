% THIS CODE IS USED FOR TIMESECTIONS ONLY. IT COMPUTES AN IMAGE OR A VIDEO
% FOR THE ATTRIBUTE ALONG THAT TIME SECTIONS AND SAVES THE
% IMAGES/VIDEOS/MAT FILES IN THE CURRENT DIRECTORY.

% ---------------------------------------------------------------------
% Load data and set up: 
% ----------------------------------------------------------------------

clc; clear; close all; 

if ismac
    seismicVolume = load('crossline.mat');
else
    seismicVolume = load('crossline.mat');
end
seismicVolume = seismicVolume.Crossline;
seismicVolume = shiftdim(seismicVolume,1); %only if dimensions not correct

[depth, iline, xline] = size(seismicVolume);

%% ---------------------------------------------------------------------
% MAIN PARAMETERS:  
% ----------------------------------------------------------------------

ts = 350; % time section to process
Sigma = [2.5  0  0;...
         0   2.5 0;...
         0   0   1]; % The covarience matrix used
dimension = 1; % dimension along which to rotate. 1: depth, 2: inline, 3: crossline
angular_resolution = 25; % how many degrees to rotate the kernel every iteration (in degrees)
sx = 5; % size of the analysis cube (spatially. That is inline x crossline)
st = 5; % number of samples in EACH TRACE .. increase for faults, decrease for channels


%% ---------------------------------------------------------------------
% MAIN LOOP:  `
% ----------------------------------------------------------------------

buffer = floor(sx/2);
buffer_depth = floor(st/2);

rotationAxis = {'depth','iline','xline'};
vals = diag(Sigma);
name = strcat('RK_ts_',num2str(ts) ,'_res_',num2str(angular_resolution),...
    '_cov_',num2str(vals(1)),'-',num2str(vals(2)), ...
    '-',num2str(vals(3)),'_size',num2str(sx),'x',...
    num2str(st),'_r_',rotationAxis{dimension});
Video1 = VideoWriter(strcat(name,'.avi'),'Uncompressed AVI');
Video2 = VideoWriter(strcat(name,'_gray.avi'),'Uncompressed AVI');
open(Video1);
open(Video2);

% tensor to save results
result_rgb = zeros(ceil(360/angular_resolution) + 1, xline, iline, 3);

% center coords of the cube: 
centerx = ceil(sx/2);
centert = ceil(st/2);
centers = [centert, centerx, centerx];

index = 0;
for theta_deg = 0:angular_resolution:180 % due to symmetry, 180 is sufficient
    index = index + 1;
    theta_rad = (2*pi/360) * theta_deg;
    R     = get3DRotationMatrix( theta_rad,dimension);
    R_inv = get3DRotationMatrix(-theta_rad,dimension);

    % EC1, EC2, and EC3:
    CC1 = zeros(iline,xline); % zeros(size(seismicVolume));
    CC2 = zeros(iline,xline); % zeros(size(seismicVolume));
    CC3 = zeros(iline,xline); % zeros(size(seismicVolume));
    
    %============================================
    % COMPUTE COHERENCE CUBE (CC)
    %============================================
    
    if ~isempty(Sigma)
        G = zeros(st,sx,sx); % Gaussian preprocessing tensor
        for ii = 1:st
            for jj = 1:sx
                for kk = 1:sx
                    G(ii,jj,kk) = exp(-0.5*([ii,jj,kk]-centers)/(R*Sigma^3*R_inv)*([ii,jj,kk]-centers)');
                    if G(ii,jj,kk) == inf
                        G(ii,jj,kk) = 0;
                    end
                end
            end
        end
    else
        G = ones(st,sx,sx);
    end
    
    % Normalize G to have a maximum value of 1:
    G = G./max(G(:));
    
    
    i = ts;
    for j = 1+buffer : iline-buffer
        
        disp(['Angle: ', num2str(index), '/', num2str(180/angular_resolution + 1), ...
            ' | Inline: ',num2str(100*(j-1-buffer)/(iline-2*buffer-1)),'%']);
        
        for k = 1+buffer : xline-buffer
            % crop analysis cube
            aCube = squeeze(seismicVolume(i-buffer_depth:i+buffer_depth, j-buffer:j+buffer, k-buffer:k+buffer));
            aCube = aCube.*G;
            
            % compute coherence along time mode
            D1 = tens2mat(aCube,1,[2 3]);
            D_tilde1 = D1 - ones(st,1)*mean(D1,1);
            C = D_tilde1'*D_tilde1;
            lambda = max(eig(C));
            CC1(j,k)=  lambda/trace(C);
            
            % compute coherence along inline mode
            D2 = tens2mat(aCube,2,[1 3]);
            D_tilde2 = D2 - ones(sx,1)*mean(D2,1);
            C = D_tilde2'*D_tilde2;
            lambda = max(eig(C));
            CC2(j,k)=  lambda/trace(C);
            
            % compute coherence along crossline mode
            D3 = tens2mat(aCube,3,[1 2]);
            D_tilde3 = D3 - ones(sx,1)*mean(D3,1);
            C = D_tilde3'*D_tilde3;
            lambda = max(eig(C));
            CC3(j,k)=  lambda/trace(C);
        end
    end
    
    % concatenate results:
    im  = cat(3, CC1, CC2, CC3);
    
    close;
    figure;
    im = rot90(im,3);
    im = fliplr(im);
    
    result_rgb(index,:,:,:) = im;
    
    
    
    image = insertText(im,[1,1],['theta: ', num2str(theta_deg)], ...
        'FontSize',20,'BoxOpacity',0.8,'TextColor','Black');
    imshow(image,[]);
    writeVideo(Video1,getframe);
    writeVideo(Video1,getframe);
    writeVideo(Video1,getframe);
    
    %     grayIm = sqrt(squeeze(im(:,:,1)).^2 + squeeze(im(:,:,2)).^2 + squeeze(im(:,:,3)).^2);
    grayIm = rgb2gray(im);
    grayIm = insertText(grayIm,[1,1],['theta: ', num2str(theta_deg)], ...
        'FontSize',20,'BoxOpacity',0.8,'TextColor','Black');
    grayIm = imNormalize(grayIm,2);
    imshow(grayIm);
    writeVideo(Video2,getframe);
    writeVideo(Video2,getframe);
    writeVideo(Video2,getframe);
    clf;
    close;
    
end

close(Video1);
close(Video2);

filename = strcat(name,'.mat');
save(filename,'result_rgb','-v7.3');


