%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gao Chenzhong
%   Beijing Key Laboratory of Fractional Signals and Systems,
%   Multi-Dimensional Signal and Information Processing Institute,
%   School of Information and Electronics, Beijing Institute of Technology
% Contact: gao-pingqi@qq.com

% Core Process of the Multisource/Multimodal Image Matching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cor1,cor2] = main_image_matching(I1, I2, parameters, ...
    int_flag, rot_flag, scl_flag, par_flag, show_flag)
%% Parameters
nOctaves1 = parameters.nOctaves1;  % Gaussian pyramid octave number, default: 3
nOctaves2 = parameters.nOctaves2;
nLayers   = parameters.nLayers;    % Gaussian pyramid layer number, default: 4
G_resize  = parameters.G_resize;   % Gaussian pyramid downsampling ratio, default: 2
G_sigma   = parameters.G_sigma;    % Gaussian blurring standard deviation, default: 1.6
key_type  = parameters.key_type;
radius    = parameters.radius;     % Local non-maximum suppression radius, default: 2
thresh    = parameters.thresh;     % Keypoints response threshold, default: 50
Npoint    = parameters.Npoint;     % Keypoints number threshold, default: 5000
patchsize = parameters.patch_size; % GGLOH patchsize, default: 72 or 96
NBA       = parameters.NBA;        % GGLOH localtion division, default: 12
NBO       = parameters.NBO;        % GGLOH orientation division, default: 12
Error     = parameters.Error;      % Outlier removal pixel loss
K         = parameters.K;          % Outlier removal repetition times

%% Keypoints detection
% LNMS Parameters of keypoints detection
% r1 = radius; r2 = radius; 
ratio = sqrt((size(I1,1)*size(I1,2))/(size(I2,1)*size(I2,2)));
if ratio>=1
    r2 = radius; r1 = round(radius*ratio);
else
    r1 = radius; r2 = round(radius/ratio);
end
tic,keypoints_1 = Detect_Keypoint(I1,6,thresh,r1,Npoint,nOctaves1,G_resize,key_type,show_flag);
    t(1)=toc; fprintf(['已完成参考图像特征点检测，用时 ',num2str(t(1)),'s\n']);
              fprintf([' Done keypoints detection of reference image, time: ',num2str(t(1)),'s\n']);
    if isempty(keypoints_1)
        cor1 = []; cor2 = []; return
    end
tic,keypoints_2 = Detect_Keypoint(I2,6,thresh,r2,Npoint,nOctaves2,G_resize,key_type,show_flag);
    t(2)=toc; fprintf(['已完成待配准图像特征点检测，用时 ',num2str(t(2)),'s\n']);
              fprintf([' Done keypoints detection of sensed image, time: ',num2str(t(2)),'s\n\n']);
    if isempty(keypoints_2)
        cor1 = []; cor2 = []; return
    end

%% Keypoints description
tic,descriptors_1 = Multiscale_Descriptor(I1,keypoints_1,patchsize,NBA,NBO,...
    nOctaves1,nLayers,G_resize,G_sigma,int_flag,rot_flag,par_flag);
    t(3)=toc; fprintf(['已完成参考图像描述符建立，用时 ',num2str(t(3)),'s\n']);
              fprintf([' Done keypoints description of reference image, time: ',num2str(t(3)),'s\n']);
tic,descriptors_2 = Multiscale_Descriptor(I2,keypoints_2,patchsize,NBA,NBO,...
    nOctaves2,nLayers,G_resize,G_sigma,int_flag,rot_flag,par_flag);
    t(4)=toc; fprintf(['已完成待配准图像描述符建立，用时 ',num2str(t(4)),'s\n']);
              fprintf([' Done keypoints description of sensed image, time: ',num2str(t(4)),'s\n\n']);

%% Keypoints matching
tic,[cor1,cor2] = Multiscale_Matching(descriptors_1,descriptors_2,...
    nOctaves1,nOctaves2,nLayers,Error,scl_flag,par_flag,K);
    t(5)=toc; fprintf(['已完成特征点匹配，用时 ',num2str(t(5)),'s\n']);
              fprintf([' Done keypoints matching, time: ',num2str(t(5)),'s\n\n']);
    if show_flag
        Show_Matches(I1,I2,cor1,cor2,1);
    end