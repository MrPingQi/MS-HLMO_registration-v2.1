%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gao Chenzhong
%   Beijing Key Laboratory of Fractional Signals and Systems,
%   Multi-Dimensional Signal and Information Processing Laboratory,
%   School of Information and Electronics, Beijing Institute of Technology
% Contact: gao-pingqi@qq.com
% 
% including:
% 1-Stage Universal Multisource/Multimodal Images Automatic Registration
% 2-Stage Accurate  Multisource/Multimodal Images Automatic Registration
% 3-Stage Large Multisource/Multimodal Images Automatic Global Registration
% 3-Stage Large Multisource/Multimodal Images Automatic Local  Registration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all; clear; clc;

%% Registration framework
% method = '1-stage';
method = '2-stage';
% method = '3-stage-1';
% method = '3-stage-2';
iter = 1; % Iteration is performable in stage-2
batchsize = 512; % If '3-stage-1or2' is selected, big data are further segmented for local registration
global_flag = 1; % If '3-stage-2' is selected, do you still want global output of the big data (RAM alert!)
save_path = '.\save_image\';
%% Is there any obvious intensity difference (multi-modal)
int_flag = 1; % yes:1, no:0
%% Is there any obvious rotation difference
rot_flag = 1;
%% Is there any obvious scale difference
scl_flag = 1;
%% Do you want parallel computing in multi-scale strategy
par_flag = 1;
%% What kind of feature point do you want as the keypoint
% key_type = 'Harris';
key_type = 'ShiTomasi';
% key_type = 'PC-Harris';
% key_type = 'PC-ShiTomasi';
%% What spatial transformation model do you need at the end
% trans_form = 'similarity';
trans_form = 'affine';
% trans_form = 'projective';
% trans_form = 'polynomial';
poly_order = 2;
%% What image pair output form do you need at the end
% out_form = 'Reference';
out_form = 'Union';
% out_form = 'Inter';
% out_form = 'Geo';
%% Do you want the resolution of sensed image be changed to the reference
chg_scale = 1;
%% Do you want the visualization of registration results
Is_flag = 1; % Visualization show
I3_flag = 1; % Overlap form
I4_flag = 1; % Mosaic form
%% Image size fitting to deal with big data and reduce RAM burden
max_size_1 = 4096; min_size_1 = 256;
max_size_2 = 512; min_size_2 = 64;

%% Algorithm parameters
parameters.nOctaves1  = 3;    % Gaussian pyram id octave number, default: 3
parameters.nOctaves2  = 3;
parameters.nLayers    = 4;    % Gaussian pyramid layer number, default: 4
parameters.G_resize   = 1.2;  % Gaussian pyramid downsampling ratio, default: 2
parameters.G_sigma    = 1.6;  % Gaussian blurring standard deviation, default: 1.6
parameters.thr        = 0;
parameters.radius     = 2;    % Local non-maximum suppression radius, default: 2
parameters.N          = 5000; % Keypoints number threhold
parameters.patch_size = 72;   % GGLOH patchsize, default: 72 or 96
parameters.NBS        = 12;   % GGLOH localtion division, default: 12
parameters.NBO        = 12;   % GGLOH orientation division, default: 12
parameters.Error      = 5;    % Outlier removal pixel loss
parameters.K          = 1;    % Rxperimental repetition times

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
parameters.poly_order = poly_order;
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

%%
% [image_1,file1,DataInfo1] = Readimage;
% [image_2,file2,DataInfo2] = Readimage;
[image_1,~,DataInfo1] = Readimage(file1);
[image_2,~,DataInfo2] = Readimage(file2);
parameters.DataInfo1 = DataInfo1;
parameters.DataInfo2 = DataInfo2;

if par_flag
    if isempty(gcp('nocreate'))
        parpool(maxNumCompThreads);
    end
end
warning off
    fprintf('\n** Registration starts, have fun\n\n');

outputs = main_process(parameters,image_1,image_2);

resample1 = outputs.resample1; cor1_r = outputs.cor1;
resample2 = outputs.resample2; cor2_r = outputs.cor2;
[I1_s,~] = Preproscessing(image_1,resample1,[]);
[I2_s,~] = Preproscessing(image_2,resample2,[]); 
matchment = Show_Matches(I1_s,I2_s,cor1_r*resample1,cor2_r*resample2,1); pause(0.01)

%% Image transformation
tic
[I1_r,I2_r,I1_rs,I2_rs,I3,I4,pos] = Transformation(image_1,image_2,cor1_r,cor2_r,...
	trans_form,out_form,chg_scale,Is_flag,I3_flag,I4_flag,poly_order);
if ~isempty(DataInfo1) && ~isempty(DataInfo1.SpatialRef)
    pos(3:4) = pos(3:4)+1;
    GeoInfo2 = Create_GeoInfo(I2_r,pos,DataInfo1);
    if strcmpi(out_form,'Geo') || strcmpi(out_form,'Reference') || strcmpi(trans_form,'polynomial')
        GeoInfo1 = DataInfo1.SpatialRef;
    else
        [rows,cols,~] = size(I1_r);
        GeoInfo1 = GeoInfo2; GeoInfo1.RasterSize = [rows,cols];
    end
else
    GeoInfo1 = []; GeoInfo2 = [];
end
str=['Done: Image tranformation, time cost: ',num2str(toc),'s\n\n']; fprintf(str);
figure,imshow(I3); title('Overlap Form'); pause(0.01)
figure,imshow(I4); title('Mosaic Form'); pause(0.01)

%% Save results
Date = datestr(now,'yyyy-mm-dd_HH-MM-SS__');
correspond = cell(2,1); correspond{1} = cor1_r; correspond{2} = cor2_r;
str=[save_path,Date,'0 correspond','.mat']; Imwrite(correspond,str,[],[]);
if exist('matchment')
    if ~isempty(matchment)
        if isvalid(matchment)
            str=[save_path,Date,'0 Matching Result','.jpg']; saveas(matchment,str);
        end
    end
end
switch out_form
    case 'Reference'
        str=[save_path,Date,'1 Reference Image','.tif'];      Imwrite(image_1,str,GeoInfo1,DataInfo1);
        if Is_flag
            [I1_s,~] = Preproscessing(image_1,1,[]); 
            str=[save_path,Date,'3 Reference Image Show','.jpg']; Imwrite(I1_s,str);
        end
    otherwise
        str=[save_path,Date,'1 Reference Image','.tif'];      Imwrite(I1_r,str,GeoInfo1,DataInfo1);
        str=[save_path,Date,'3 Reference Image Show','.jpg']; Imwrite(I1_rs,str);
end
str=[save_path,Date,'2 Registered Image','.tif'];             Imwrite(I2_r,str,GeoInfo2,DataInfo1);
str=[save_path,Date,'4 Registered Image Show','.jpg'];        Imwrite(I2_rs,str);
str=[save_path,Date,'5 Overlap of results','.jpg'];           Imwrite(I3,str);
str=[save_path,Date,'6 Mosaic of results','.jpg'];            Imwrite(I4,str);
    str='The results are saved in the save_image folder.\n\n'; fprintf(str);
