%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gao Chenzhong
%   Beijing Key Laboratory of Fractional Signals and Systems,
%   Multi-Dimensional Signal and Information Processing Institute,
%   School of Information and Electronics, Beijing Institute of Technology
% Contact: gao-pingqi@qq.com
% 
% including:
% 1-Stage Universal Multisource/Multimodal Images Automatic Registration
% 2-Stage Accurate  Multisource/Multimodal Images Automatic Registration
% 3-Stage Large Multisource/Multimodal Images Automatic Global Registration
% 3-Stage Large Multisource/Multimodal Images Automatic Local  Registration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc; clear;
addpath('functions','func_Math','func_Geo')

%% Registration framework
% method = '1-stage';
% method = '2-stage';
% method = '3-stage-1';
method = '3-stage-2';
iter = 1;  % Iteration is performable in stage-2
batchsize = 1024;  % If '3-stage-1or2' is selected, big data are further segmented for local registration
global_flag = 1;  % If '3-stage-2' is selected, do you still want global output of the big data (RAM alert!)
save_path = '.\save_image\';
%% Is there any obvious intensity difference (multi-modal)
int_flag = 1;  % yes:1, no:0
%% Is there any obvious rotation difference
rot_flag = 0;
%% Is there any obvious scale difference
scl_flag = 1;
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
%% What image pair output form do you need at the end
% out_form = 'Reference';
% out_form = 'Union';
% out_form = 'Inter';
out_form = 'Geo';
%% Do you want the resolution of sensed image be changed to the reference
chg_scale = 1;
%% Do you want the visualization of registration results
Is_flag = 1;  % Visualization show
I3_flag = 1;  % Overlap form
I4_flag = 1;  % Mosaic form
%% Image size fitting to deal with big data and reduce RAM burden
max_size_1 = 4096; min_size_1 = 256;
max_size_2 = 1024; min_size_2 = 64;

%% Algorithm parameters
parameters.nOctaves1  = 3;    % Gaussian pyramid octave number, default: 3
parameters.nOctaves2  = 3;
parameters.nLayers    = 4;    % Gaussian pyramid layer number, default: 4
parameters.G_resize   = 2;    % Gaussian pyramid downsampling ratio, default: 2
parameters.G_sigma    = 1.6;  % Gaussian blurring standard deviation, default: 1.6
parameters.radius     = 2;    % Local non-maximum suppression radius, default: 2
parameters.thresh     = 0;    % Keypoints response threshold, default: 0
parameters.Npoint     = 5000; % Keypoints number threshold
parameters.patch_size = 72;   % GGLOH patchsize, default: 72 or 96
parameters.NBA        = 12;   % GGLOH localtion division, default: 12
parameters.NBO        = 12;   % GGLOH orientation division, default: 12
parameters.Error      = 5;    % Outlier removal pixel loss
parameters.K          = 1;    % Outlier removal repetition times
parameters.disp       = 1;    % Display intermediate results

%%
parameters.method = method;
parameters.iter = iter;
parameters.batchsize = batchsize;
parameters.global_flag = global_flag;
parameters.int_flag = int_flag;
parameters.rot_flag = rot_flag;
parameters.scl_flag = scl_flag;
parameters.par_flag = par_flag;
parameters.key_type = key_type;
parameters.trans_form = trans_form;
parameters.out_form = out_form;
parameters.chg_scale = chg_scale;
parameters.Is_flag = Is_flag;
parameters.I3_flag = I3_flag;
parameters.I4_flag = I4_flag;
parameters.max_size_1 = max_size_1;
parameters.min_size_1 = min_size_1;
parameters.max_size_2 = max_size_2;
parameters.min_size_2 = min_size_2;
parameters.save_path = save_path;
addpath(genpath(pwd));
if (exist(save_path,'dir')==0) % If file folder does not exist
    mkdir(save_path);
end

%% Registration process
[image_1,file1,DataInfo1] = Readimage;
[image_2,file2,DataInfo2] = Readimage;
% [image_1,~,DataInfo1] = Readimage(file1);
% [image_2,~,DataInfo2] = Readimage(file2);
parameters.DataInfo1 = DataInfo1;
parameters.DataInfo2 = DataInfo2;

if par_flag && isempty(gcp('nocreate'))
    parpool(maxNumCompThreads);  % Start parallel computing, time needed
end
fprintf('\n*开始图像配准，请耐心等待...\n Image registration starts, please be patient...\n\n');
warning off; t = []; tic

outputs = main_process(parameters,image_1,image_2);  % Main process

t(1)=toc;
resample1 = outputs.resample1; cor1 = outputs.cor1;
resample2 = outputs.resample2; cor2 = outputs.cor2;
[I1_s,~] = Preproscessing(image_1,resample1,[]);
[I2_s,~] = Preproscessing(image_2,resample2,[]);
matchment = Show_Matches(I1_s,I2_s,cor1*resample1,cor2*resample2,1);

%% Image transformation
tic,[I1_r,I2_r,I1_rs,I2_rs,I3,I4,~,pos] = Transformation(image_1,image_2,...
    cor1,cor2,trans_form,out_form,chg_scale,Is_flag,I3_flag,I4_flag);
    if ~isempty(DataInfo1) && ~isempty(DataInfo1.SpatialRef)
        pos(3:4) = pos(3:4)+1;
        GeoInfo2 = Create_GeoInfo(I2_r,pos,DataInfo1);
        if strcmpi(out_form,'Geo') || strcmpi(out_form,'Reference')
            GeoInfo1 = DataInfo1.SpatialRef;
        else
            [rows,cols,~] = size(I1_r);
            GeoInfo1 = GeoInfo2; GeoInfo1.RasterSize = [rows,cols];
        end
    else
        GeoInfo1 = []; GeoInfo2 = [];
    end
    t(2)=toc; fprintf(['已完成图像变换，用时 ',num2str(t(2)),'s\n']);
              fprintf([' Done image transformation，time: ',num2str(t(2)),'s\n\n']);
    figure,imshow(I3),title('Overlap Form'); drawnow
    figure,imshow(I4),title('Mosaic Form'); drawnow

%%
T=num2str(sum(t)); fprintf(['*已完成图像配准，总用时 ',T,'s\n']);
                   fprintf([' Done image registration, total time: ',T,'s\n\n']);

%% Save results
Date = datestr(now,'yyyy-mm-dd_HH-MM-SS__'); tic
cors = {cor1;cor2}; Imwrite(cors,[save_path,Date,'0 corresponds','.mat']);
if exist('matchment','var') && ~isempty(matchment) && isvalid(matchment)
    saveas(matchment,[save_path,Date,'0 Matching result','.jpg']);
end
if strcmpi(out_form,'reference')
    Imwrite(image_1,[save_path,Date,'1 Reference image','.tif'],GeoInfo1,DataInfo1);
    if Is_flag
        [I1_s,~] = Preproscessing(image_1,1,[]); 
        Imwrite(I1_s,[save_path,Date,'3 Reference image show','.jpg']);
    end
else
    Imwrite(I1_r ,[save_path,Date,'1 Reference image','.tif'],GeoInfo1,DataInfo1);
    Imwrite(I1_rs,[save_path,Date,'3 Reference image show','.jpg']);
end
Imwrite(I2_r ,[save_path,Date,'2 Registered image','.tif'],GeoInfo2,DataInfo1);
Imwrite(I2_rs,[save_path,Date,'4 Registered image show','.jpg']);
Imwrite(I3   ,[save_path,Date,'5 Overlap of results','.jpg']);
Imwrite(I4   ,[save_path,Date,'6 Mosaic of results','.jpg']);
T=num2str(toc); fprintf(['配准结果已经保存在程序根目录下的save_image文件夹中，用时',T,'s\n']);
                fprintf([' Registration results are saved in the save_image folder, time: ',T,'s\n']);