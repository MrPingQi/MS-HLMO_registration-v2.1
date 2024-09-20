%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gao Chenzhong
%   Beijing Key Laboratory of Fractional Signals and Systems,
%   Multi-Dimensional Signal and Information Processing Institute,
%   School of Information and Electronics, Beijing Institute of Technology
% Contact: gao-pingqi@qq.com

% Main Process of Multisource/Multimodal Images Automatic Registration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outputs = main_process(parameters,image_1,image_2)
if exist('parameters','var') && ~isempty(parameters)
    func_flag = 1;
    method     = parameters.method;
    iter       = parameters.iter;
    batchsize  = parameters.batchsize;
    global_flag = parameters.global_flag;
    int_flag   = parameters.int_flag;
    rot_flag   = parameters.rot_flag;
    scl_flag   = parameters.scl_flag;
    par_flag   = parameters.par_flag;
    key_type   = parameters.key_type;
    trans_form = parameters.trans_form;
    Is_flag    = parameters.Is_flag;
    max_size_1 = parameters.max_size_1;
    min_size_1 = parameters.min_size_1;
    max_size_2 = parameters.max_size_2;
    min_size_2 = parameters.min_size_2;
    DataInfo1  = parameters.DataInfo1;
    DataInfo2  = parameters.DataInfo2;
    save_path  = parameters.save_path;
else
    close all; clear; clc;
    func_flag = 0;
    method = '2-stage';
    iter = 1;
    batchsize = 1024;
    global_flag = 1;
    int_flag = 1;
    rot_flag = 1;
    scl_flag = 1;
    par_flag = 1;
    key_type = 'ShiTomasi';  % Harris, ShiTomasi, PC-Harris, PC-ShiTomasi
    trans_form = 'affine';  % similarity, affine, projective
    Is_flag = 1;  % Visualization show
    max_size_1 = 4096;
    min_size_1 = 256;
    max_size_2 = 1024;
    min_size_2 = 64;
    
    parameters.nOctaves1  = 3;    % Gaussian pyramid octave number, default: 3
    parameters.nOctaves2  = 3;
    parameters.nLayers    = 4;    % Gaussian pyramid layer number, default: 4
    parameters.G_resize   = 2;    % Gaussian pyramid downsampling ratio, default: 2
    parameters.G_sigma    = 1.6;  % Gaussian blurring standard deviation, default: 1.6
    parameters.thresh     = 0;    % Keypoints response threshold
    parameters.radius     = 2;    % Local non-maximum suppression radius, default: 2
    parameters.Npoint     = 5000; % Keypoints number threshold
    parameters.patch_size = 72;  % GGLOH patchsize, default: 72 or 96
    parameters.NBA        = 12;   % GGLOH localtion division, default: 12
    parameters.NBO        = 12;   % GGLOH orientation division, default: 12
    parameters.Error      = 5;    % Outlier removal pixel loss
    parameters.K          = 1;    % Outlier removal repetition times
    parameters.disp       = 1;    % Display intermediate results
    
    [image_1,~,DataInfo1] = Readimage;
    [image_2,~,DataInfo2] = Readimage;
end

%% Initialization
trans1_t = cell(8,2); trans2_t = cell(8,2);
trans1_t(1:end,1) = {[1,0,0;0,1,0;0,0,1]}; trans1_t(1:end,2) = {'affine'};
trans2_t(1:end,1) = {[1,0,0;0,1,0;0,0,1]}; trans2_t(1:end,2) = {'affine'};

%% Image preproscessing
[M1,N1,~] = size(image_1); [M2,N2,B2] = size(image_2);
[~,resample1] = Deal_Extreme(image_1,min_size_1,max_size_1,0);
[I1_s,I1_o] = Preproscessing(image_1,resample1);
figure,imshow(I1_s),title('Reference image base'); drawnow; clear I1_s
if ~strcmpi(method,'3-stage-2')
    clear image_1
end
[~,resample2] = Deal_Extreme(image_2,min_size_1,max_size_1,0);
[I2_s,I2_o] = Preproscessing(image_2,resample2);
figure,imshow(I2_s),title('Sensed Image base'); drawnow; clear I2_s
if ~strcmpi(method,'3-stage-2')
    clear image_2
end
outputs.resample1 = resample1; outputs.resample2 = resample2;
trans1_t{1,1} = [resample1,0,0; 0,resample1,0; 0,0,1];
trans2_t{1,1} = [resample2,0,0; 0,resample2,0; 0,0,1];

I1 = I1_o; I2 = I2_o;
%%
% if ~isempty(DataInfo1) && ~isempty(DataInfo1.SpatialRef) &&...
%    ~isempty(DataInfo2) && ~isempty(DataInfo2.SpatialRef)
%     Cs = [1,1; N2,1; 1,M2; N2,M2];
%     % Sensed image pixel to world
%     points = Pix2world(Cs, DataInfo2);
%     % World to reference image pixel
%     Cs = World2pix(points, DataInfo1);
%     A_x = max(1,min(Cs(:,1))); A_y = max(1,min(Cs(:,2)));
%     C_x = min(N1,max(Cs(:,1))); C_y = min(M1,max(Cs(:,2)));
%     if C_x-A_x>min_size_1 && C_y-A_y>min_size_1
%         tiepoints1 = Cs; clear Cs A_x A_y C_x C_y;
%         [M,N,~] = size(I2_o);
%         cor2_t = [1,1; N,1; 1,M; N,M]; cor1_t = tiepoints1*resample1;
%         [I1,I2,trans_1,trans_2] = ...
%             Transform_temp(I1_o,I2_o,cor1_t,cor2_t,'projective');
%         trans1_t(2,:) = trans_1;
%         trans2_t(2,:) = trans_2;
%         rot_flag = 0; scl_flag = 0;
%     else
%         fprintf('No overlap area or Geo infomation incorrect');
%         I1 = I1_o; I2 = I2_o;
%     end
% else
%     I1 = I1_o; I2 = I2_o;
% end
% 
% if manu_flag
%     [~,I1,roi1] = Select_ROI(Visual(I1),I1,1);
%     [~,I2,roi2] = Select_ROI(Visual(I2),I2,1);
%     trans1_t{3,1} = [1,0,0; 0,1,0; -roi1(1)+1,-roi1(2)+1,1];
%     trans2_t{3,1} = [1,0,0; 0,1,0; -roi2(1)+1,-roi2(2)+1,1];
%     
%     [trans_t,~,~] = Manual_Register(Visual(I1),Visual(I2),'projective');
%     trans1_t{4,1} = [1,0,0;0,1,0;0,0,1];
%     trans2_t(4,:) = trans_t;
%     t_form = projective2d(trans2_t{4,1});
%     I2 = imwarp(I2,t_form,'OutputView',imref2d(size(I1)));
% end

[I1,resample1_t] = Deal_Extreme(I1,min_size_2,max_size_2,1);
[I2,resample2_t] = Deal_Extreme(I2,min_size_2,max_size_2,1);
trans1_t{5,1} = [resample1_t,0,0; 0,resample1_t,0; 0,0,1];
trans2_t{5,1} = [resample2_t,0,0; 0,resample2_t,0; 0,0,1];

%%
if par_flag && isempty(gcp('nocreate'))
    parpool(maxNumCompThreads);  % Start parallel computing, time needed
end
if ~func_flag
    fprintf('\n*开始图像配准，请耐心等待...\n Image registration starts, please be patient...\n\n');
    warning off;
end


%% Stage 1
fprintf('* Stage 1\n');
[cor1,cor2] = main_image_matching(I1,I2,parameters,...
    int_flag,rot_flag,scl_flag,par_flag,0);

[cor1_1,~] = XY_Transform(cor1,trans1_t(2:end,:),-1);
[cor2_1,~] = XY_Transform(cor2,trans2_t(2:end,:),-1);

[I1,I2,trans_t1,trans_t2] = ...
    Transform_temp(I1_o,I2_o,cor1_1,cor2_1,'affine');
figure,imshow(Visual(I1)/2+Visual(I2)/2,[]);
trans1_t(6,:) = trans_t1; cor1_o = cor1_1;
trans2_t(6,:) = trans_t2; cor2_o = cor2_1;


%% Stage 2
if ~strcmpi(method,'1-stage') 
fprintf('* Stage 2\n');
for k = 1:iter
if iter>1
    fprintf(['Iteration = ',num2str(k),'\n']);
end
[I1,resample1_t] = Deal_Extreme(I1,min_size_2,max_size_2,1);
[I2,resample2_t] = Deal_Extreme(I2,min_size_2,max_size_2,1);
trans1_t{7,1} = [resample1_t,0,0; 0,resample1_t,0; 0,0,1];
trans2_t{7,1} = [resample2_t,0,0; 0,resample2_t,0; 0,0,1];

[cor1,cor2] = main_image_matching(I1,I2,parameters,...
    int_flag,0,0,par_flag,0);

[cor1,~] = XY_Transform(cor1,trans1_t(7,:),-1);
[cor2,~] = XY_Transform(cor2,trans2_t(7,:),-1);
[cor1_2,~] = XY_Transform(cor1,trans_t1,-1);
[cor2_2,~] = XY_Transform(cor2,trans_t2,-1);

[I1,I2,trans_t1,trans_t2] = ...
    Transform_temp(I1_o,I2_o,cor1_2,cor2_2,trans_form);
figure,imshow(Visual(I1)/2+Visual(I2)/2,[]);
end
trans1_t(8,:) = trans_t1; cor1_o = cor1_2;
trans2_t(8,:) = trans_t2; cor2_o = cor2_2;
end


%% Stage 3.1
if strcmpi(method,'3-stage-1')
fprintf('* Stage 3\n');
expand = ceil(batchsize/4);
iter = 1;
for k = 1:iter
if iter>1
    fprintf(['Iteration = ',num2str(k),'\n']);
end
[M,N] = size(I1);
batch_x = max(round(N/batchsize),1); XX = [((0:batch_x-1)*batchsize),N];
batch_y = max(round(M/batchsize),1); YY = [((0:batch_y-1)*batchsize),M];
cor1_t = []; cor2_t = []; nbatch = num2str(batch_x*batch_y); n = 0;
for i=1:batch_y
    yy1 = YY(i)+1; y1 = max(yy1-expand,1);
    yy2 = YY(i+1); y2 = min(yy2+expand,M);
    for j=1:batch_x
        xx1 = XX(j)+1; x1 = max(xx1-expand,1);
        xx2 = XX(j+1); x2 = min(xx2+expand,N);
        n = n+1; fprintf([' ',num2str(n),'/',nbatch,'\n']);
        trans_cut1 = {[1,0,0;0,1,0;-x1+1,-y1+1,1],'affine'};
        trans_cut2 = {[1,0,0;0,1,0;-x1+1,-y1+1,1],'affine'};
        [cor1,cor2] = main_image_matching(I1(y1:y2,x1:x2),....
            I2(y1:y2,x1:x2),parameters,int_flag,0,0,par_flag,0);
        [cor1,~] = XY_Transform(cor1,trans_cut1,-1);
        [cor2,~] = XY_Transform(cor2,trans_cut2,-1);
        if size(cor1,1)>20
            cor1_t = [cor1_t; cor1];
            cor2_t = [cor2_t; cor2];
        end
    end
end
matches = [cor1_t,cor2_t];
[~,index1,~] = unique(matches(:,1:2),'rows');
matches = matches(index1,:);
[~,index2,~] = unique(matches(:,3:4),'rows');
matches = matches(index2,:);
cor1 = matches(:,1:2); cor2 = matches(:,3:4);

[cor1_3,~] = XY_Transform(cor1,trans_t1,-1);
[cor2_3,~] = XY_Transform(cor2,trans_t2,-1);

[I1,I2,trans_t1,trans_t2] = ...
    Transform_temp(I1_o,I2_o,cor1_3,cor2_3,trans_form);
figure,imshow(Visual(I1)/2+Visual(I2)/2,[]);
end
cor1_o = cor1_3; cor2_o = cor2_3;
end


%% Stage 3.2
if strcmpi(method,'3-stage-2')
fprintf('* Stage 3\n');
expand = ceil(batchsize/4);
[M,N] = size(I1);
% XX,YY 用以保存大图像列、行切割界限
batch_x = max(round(N/batchsize),1); XX = [((0:batch_x-1)*batchsize),N];
batch_y = max(round(M/batchsize),1); YY = [((0:batch_y-1)*batchsize),M];

dX = -trans1_t{8,1}(3,1); dY = -trans1_t{8,1}(3,2);
% XXr,YYr 用以保存大图像原始数据对应切割界限
XXr = round((XX+dX)/resample1); XXr = min(XXr,N1);
YYr = round((YY+dY)/resample1); YYr = min(YYr,M1);

if global_flag
    I2_r = zeros(M1,N1,B2);
end
cor1_3 = []; cor2_3 = []; nbatch = num2str(batch_x*batch_y); n = 0;
for i=1:batch_y
    yy1 = YY(i)+1; y1 = max(yy1-expand,1);
    yy2 = YY(i+1); y2 = min(yy2+expand,M);
    for j=1:batch_x
        xx1 = XX(j)+1; x1 = max(xx1-expand,1);
        xx2 = XX(j+1); x2 = min(xx2+expand,N);
        n = n+1; fprintf([' ',num2str(n),'/',nbatch,'\n']);
        trans_cut1 = {[1,0,0;0,1,0;-x1+1,-y1+1,1],'affine'};
        trans_cut2 = {[1,0,0;0,1,0;-x1+1,-y1+1,1],'affine'};
        [cor1,cor2] = main_image_matching(I1(y1:y2,x1:x2),....
            I2(y1:y2,x1:x2),parameters,int_flag,0,0,par_flag,0);
        if size(cor1,1)>20
            [cor2,~] = XY_Transform(cor2,trans_cut2,-1);
            [cor2,~] = XY_Transform(cor2,trans2_t(8,:),-1);
            cor2_3 = [cor2_3; cor2];
            [cor2_t,~] = XY_Transform(cor2,trans2_t(1,:),-1);
            
            cor1_t = [cor1(:,1)-(xx1-x1),cor1(:,2)-(yy1-y1)];
            [cor1_t,~] = XY_Transform(cor1_t,trans1_t(1,:),-1);  % 注意image_1取切块方式
            [cor1,~] = XY_Transform(cor1,trans_cut1,-1);
            [cor1,~] = XY_Transform(cor1,trans1_t(8,:),-1);
            cor1_3 = [cor1_3; cor1];
            
            Show_Matches(I1_o,I2_o,cor1,cor2,1);
            
            xxr1 = XXr(j)+1; xxr2 = XXr(j+1);
            yyr1 = YYr(i)+1; yyr2 = YYr(i+1);
            
            I1_r_t = image_1(yyr1:yyr2,xxr1:xxr2,:);
            [~,I2_r_t,I1_rs_t,I2_rs_t,~,~,~,~] = Transformation(I1_r_t,image_2,...
                cor1_t,cor2_t,trans_form,'Reference',1,Is_flag,0,0);
            if global_flag
                I2_r(yyr1:yyr2,xxr1:xxr2,:) = I2_r_t;
            end
            
            pos = [xxr1,yyr1,xxr2+1,yyr2+1];
            GeoInfo = Create_GeoInfo(I2_r_t,pos,DataInfo1);
            
            Date = datestr(now,'yyyy-mm-dd_HH-MM-SS__');
            Imwrite(I1_r_t ,[save_path,Date,'1 Reference image' ,'.tif'],GeoInfo,DataInfo1);
            Imwrite(I2_r_t ,[save_path,Date,'2 Registered image','.tif'],GeoInfo,DataInfo1);
            Imwrite(I1_rs_t,[save_path,Date,'1 Reference image' ,'.jpg']);
            Imwrite(I2_rs_t,[save_path,Date,'2 Registered image','.jpg']);
        end
    end
end
matches = [cor1_3,cor2_3];
[~,index1,~] = unique(matches(:,1:2),'rows');
matches = matches(index1,:);
[~,index2,~] = unique(matches(:,3:4),'rows');
matches = matches(index2,:);
cor1_3 = matches(:,1:2); cor2_3 = matches(:,3:4);
cor1_o = cor1_3; cor2_o = cor2_3;
end

if global_flag
    if ~isempty(DataInfo1) && ~isempty(DataInfo1.SpatialRef)
        GeoInfo = DataInfo1.SpatialRef;
    else
        GeoInfo = [];
    end
    Date = datestr(now,'yyyy-mm-dd_HH-MM-SS__');
    Imwrite(I2_r,[save_path,Date,'Registered image','.tif'],GeoInfo,DataInfo1);
end


%% Postprocessing
[cor1_r,~] = XY_Transform(cor1_o,trans1_t(1,:),-1);
[cor2_r,~] = XY_Transform(cor2_o,trans2_t(1,:),-1);
outputs.cor1 = cor1_r; outputs.cor2 = cor2_r;