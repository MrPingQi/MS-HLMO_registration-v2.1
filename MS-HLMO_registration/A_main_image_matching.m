%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gao Chenzhong
%   Beijing Key Laboratory of Fractional Signals and Systems,
%   Multi-Dimensional Signal and Information Processing Laboratory,
%   School of Information and Electronics, Beijing Institute of Technology
% Contact: gao-pingqi@qq.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc;

%% Is there any obvious intensity difference (multi-modal)
int_flag = 1; % yes:1, no:0
%% Is there any obvious rotation difference
rot_flag = 1;
%% Is there any obvious scale difference
scl_flag = 0;
%% Do you want parallel computing in multi-scale strategy
par_flag = 1;
%% What kind of feature point do you want as the keypoint
% key_type = 'Harris';
% key_type = 'ShiTomasi';
% key_type = 'PC-Harris';
key_type = 'PC-ShiTomasi';
%% What spatial transformation model do you need at the end
% trans_form = 'similarity';
trans_form = 'affine';
% trans_form = 'projective';
% trans_form = 'polynomial-n'; % n = 2,3,...
%% What transformation model do you need at the end
% output_form = 'Reference';
output_form = 'Union';
% output_form = 'Inter';

%% Parameters
G_resize = 2;  % Gaussian pyramid downsampling ratio, default: 2
nOctaves1 = 3; % Gaussian pyramid octave number, default: 3
nOctaves2 = 3; 
G_sigma = 1.6; % Gaussian blurring standard deviation, default: 1.6
nLayers = 4;   % Gaussian pyramid layer number, default: 4
radius = 2;    % Local non-maximum suppression radius, default: 2
N = 2000;      % Keypoints number threhold
patch_size = 72; % GGLOH patchsize, default: 72 or 96
NBS = 12;        % GGLOH localtion division, default: 12
NBO = 12;        % GGLOH orientation division, default: 12
Error = 5;       % Outlier removal pixel loss
K = 1;           % Rxperimental repetition times

%% Read Images
[image_1, image_2] = Readimage;

%% Image Preproscessing
resample1 = 1/1; resample2 = 1/1;
[I1_s,I1] = Preproscessing(image_1,resample1,[]);
[I2_s,I2] = Preproscessing(image_2,resample2,[]);
figure,imshow(I1_s),title('Reference image'); pause(0.01)
figure,imshow(I2_s),title('Sensed Image'); pause(0.01)
% figure,imshow(I1),title('Reference image'); pause(0.01)
% figure,imshow(I2),title('Sensed Image'); pause(0.01)

if size(resample1,2)==1
    resample1 = [resample1,resample1];
end
if size(resample2,2)==1
    resample2 = [resample2,resample2];
end

%%
warning off
    fprintf('\n** Registration starts, have fun\n\n'); ts=cputime;

%% Keypoints Detection
% LNMS Parameters of keypoints detection
% r1 = radius; r2 = radius; 
ratio = sqrt((size(I1,1)*size(I1,2))/(size(I2,1)*size(I2,2)));
if ratio>=1
    r2 = radius; r1 = round(radius*ratio);
else
    r1 = radius; r2 = round(radius/ratio);
end
tic,keypoints_1 = Detect_Keypoint(I1,6,r1,N,nOctaves1,G_resize,key_type);
    str=['Done: Keypoints detection of reference image, time cost: ',num2str(toc),'s\n']; fprintf(str);
    figure,imshow(I1_s); hold on; plot(keypoints_1(:,1)/resample1(1),keypoints_1(:,2)/resample1(2),'r+'); pause(0.01)
tic,keypoints_2 = Detect_Keypoint(I2,6,r2,N,nOctaves2,G_resize,key_type);
    str=['Done: Keypoints detection of sensed image, time cost: ',num2str(toc),'s\n\n']; fprintf(str);
    figure,imshow(I2_s); hold on; plot(keypoints_2(:,1)/resample2(1),keypoints_2(:,2)/resample2(2),'r+'); pause(0.01)
%     figure,imshow(I1); hold on; plot(keypoints_1(:,1),keypoints_1(:,2),'r+'); pause(0.01)
%     figure,imshow(I2); hold on; plot(keypoints_1(:,1),keypoints_1(:,2),'r+'); pause(0.01)

%% Keypoints Description
tic,descriptors_1 = Multiscale_Descriptor(I1,keypoints_1,patch_size,NBS,NBO,...
    nOctaves1,nLayers,G_resize,G_sigma,int_flag,rot_flag,par_flag);
    str=['Done: Keypoints description of reference image, time cost: ',num2str(toc),'s\n']; fprintf(str);
tic,descriptors_2 = Multiscale_Descriptor(I2,keypoints_2,patch_size,NBS,NBO,...
    nOctaves2,nLayers,G_resize,G_sigma,int_flag,rot_flag,par_flag);
    str=['Done: Keypoints description of sensed image, time cost: ',num2str(toc),'s\n\n']; fprintf(str);

%% Keypoints Matching
tic
if K==1
    [cor1_o,cor2_o] = Multiscale_Matching(descriptors_1,descriptors_2,...
        nOctaves1,nOctaves2,nLayers,scl_flag,par_flag);
    [cor1,cor2] = Outlier_Removal(cor1_o,cor2_o,Error);
else
    aaa = []; correspond_1 = cell(K,1); correspond_2 = cell(K,1);
    for k = 1:K
        [cor1_o,cor2_o] = Multiscale_Matching(descriptors_1,descriptors_2,...
            nOctaves1,nOctaves2,nLayers,scl_flag,par_flag);
        [cor1,cor2] = Outlier_Removal(cor1_o,cor2_o,Error);
        correspond_1{k} = cor1; correspond_2{k} = cor2; aaa = [aaa,size(cor1,1)];
    end
    [~,index] = max(aaa);
    cor1 = correspond_1{index}; cor2 = correspond_2{index};
end
    str = ['Done: Keypoints matching, time cost: ',num2str(toc),'s\n\n']; fprintf(str);

matchment = Show_Matches(I1_s,I2_s,cor1,cor2,1); pause(0.01)

%% Image transformation
tic
switch output_form
    case 'Reference'
        [I2_r,I2_rs,I3,I4] = Transform_ref...
            (image_1,image_2,cor1,cor2,trans_form);
    case 'Union'
        [I1_r,I2_r,I1_rs,I2_rs,I3,I4] = Transform_union...
            (image_1,image_2,cor1,cor2,trans_form);
    case 'Inter'
        [I1_r,I2_r,I1_rs,I2_rs,I3,I4] = Transform_inter...
            (image_1,image_2,cor1,cor2,trans_form);
end
    str=['Done: Image tranformation, time cost: ',num2str(toc),'s\n\n']; fprintf(str); tic
figure; imshow(I3,[]); title('Fusion Form'); pause(0.01)
figure; imshow(I4,[]); title('Mosaic Form'); pause(0.01)

%% Save results
if (exist('save_image','dir')==0) % If file folder does not exist
    mkdir('save_image');
end
Date = datestr(now,'yyyy-mm-dd_HH-MM-SS__');
correspond = cell(2,1); correspond{1} = cor1; correspond{2} = cor2;
str=['.\save_image\',Date,'0 correspond','.mat']; save(str,'correspond')
if isvalid(matchment)
    str=['.\save_image\',Date,'0 Matching Result','.jpg']; saveas(matchment,str);
end
switch output_form
    case 'Reference'
        str=['.\save_image\',Date,'1 Reference Image','.mat']; save(str,'image_1');
        str=['.\save_image\',Date,'1 Reference Image','.jpg']; imwrite(I1_s,str);
    otherwise
        str=['.\save_image\',Date,'1 Reference Image','.mat']; save(str,'I1_r');
        str=['.\save_image\',Date,'1 Reference Image','.jpg']; imwrite(I1_rs,str);
end
str=['.\save_image\',Date,'2 Registered Image','.mat']; save(str,'I2_r');
str=['.\save_image\',Date,'2 Registered Image','.jpg']; imwrite(I2_rs,str);
str=['.\save_image\',Date,'3 Fusion of results','.jpg']; imwrite(I3,str);
str=['.\save_image\',Date,'4 Mosaic of results','.jpg']; imwrite(I4,str);
    str='The results are saved in the save_image folder.\n\n'; fprintf(str);