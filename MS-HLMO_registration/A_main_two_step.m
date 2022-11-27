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
%% What image pair output form do you need at the end
% output_form = 'Reference';
output_form = 'Union';
% output_form = 'Inter';

%% Read Images
[image_1, image_2] = Readimage;

%% Image Preproscessing
resample1 = 1/1; resample2 = [1,1]*1/1;
[I1_s,I1_o] = Preproscessing(image_1,resample1,[]);
[I2_s,I2_o] = Preproscessing(image_2,resample2,[]);
figure,imshow(I1_s),title('Reference image'); pause(0.01)
figure,imshow(I2_s),title('Sensed Image'); pause(0.01)
% figure,imshow(I1_o),title('Reference image'); pause(0.01)
% figure,imshow(I2_o),title('Sensed Image'); pause(0.01)

if size(resample1,2)==1
    resample1 = [resample1,resample1];
end
if size(resample2,2)==1
    resample2 = [resample2,resample2];
end

%%
warning off
    fprintf('\n** Registration starts, have fun\n\n'); ts=cputime;

%% Stage 1
[cor1,cor2] = main_image_matching(I1_o,I2_o,...
    int_flag,rot_flag,scl_flag,par_flag,key_type);
Show_Matches(I1_o,I2_o,cor1,cor2,1);pause(0.01)

t_form = fitgeotrans(cor2(:,1:2),cor1(:,1:2),'similarity');
solution_t = t_form.T;
[I1,I2,dX,dY] = Transform_temp(I1_o,I2_o,solution_t);
figure,imshow(I1/2+I2/2,[]);

%% Stage 2
iter = 1;
for k = 1:iter
[cor1,cor2] = main_image_matching(I1,I2,...
    int_flag,0,0,par_flag,key_type);
t_form = fitgeotrans(cor2(:,1:2),cor1(:,1:2),trans_form);

if k~=iter
    Show_Matches(I1,I2,cor1,cor2,1);pause(0.01)
    solution_t = solution_t...
    *[1,0,0;0,1,0;-dX,-dY,1]*(t_form.T)*[1,0,0;0,1,0;dX,dY,1];
    [I1,I2,dX,dY] = Transform_temp(I1,I2,solution_t);
    figure,imshow(I1/2+I2/2,[]);
else
    cor1(:,1) = (cor1(:,1)+dX)/resample1(2); cor1(:,2) = (cor1(:,2)+dY)/resample1(1);
    cor2(:,1) = cor2(:,1)+dX; cor2(:,2) = cor2(:,2)+dY;
    cor2 = [cor2(:,1:2),ones(size(cor2,1),1)]*inv(solution_t);
    cor2(:,1) = cor2(:,1)/resample2(2); cor2(:,2) = cor2(:,2)/resample2(1);
    matchment = Show_Matches(I1_s,I2_s,cor1,cor2,1);pause(0.01)
    solution_t = solution_t...
    *[1,0,0;0,1,0;-dX,-dY,1]*(t_form.T)*[1,0,0;0,1,0;dX,dY,1];
end
end

%% Postprocessing
solution = ...
    [resample2(2),0,0; 0,resample2(1),0; 0,0,1]...
    *solution_t...
    *[1/resample1(2),0,0; 0,1/resample1(1),0; 0,0,1];

%% Image transformation
tic
switch output_form
    case 'Reference'
        [I2_r,I2_rs,I3,I4] = Transform_ref_1...
            (image_1,image_2,solution);
    case 'Union'
        [I1_r,I2_r,I1_rs,I2_rs,I3,I4] = Transform_union_1...
            (image_1,image_2,solution);
    case 'Inter'
        [I1_r,I2_r,I1_rs,I2_rs,I3,I4] = Transform_inter_1...
            (image_1,image_2,solution);
end
    str=['Done: Image tranformation, time cost: ',num2str(toc),'s\n\n']; fprintf(str); tic
figure,imshow(I3,[]); title('Fusion Form'); pause(0.01)
figure,imshow(I4,[]); title('Mosaic Form'); pause(0.01)

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
        str=['.\save_image\',Date,'1 Reference Image','.mat']; save(str,'image_1','-v7.3');
        str=['.\save_image\',Date,'3 Reference Image','.jpg']; imwrite(I1_s,str);
    otherwise
        str=['.\save_image\',Date,'1 Reference Image','.mat']; save(str,'I1_r','-v7.3');
        str=['.\save_image\',Date,'3 Reference Image','.jpg']; imwrite(I1_rs,str);
end
str=['.\save_image\',Date,'2 Registered Image','.mat']; save(str,'I2_r','-v7.3');
str=['.\save_image\',Date,'4 Registered Image','.jpg']; imwrite(I2_rs,str);
str=['.\save_image\',Date,'5 Fusion of results','.jpg']; imwrite(I3,str);
str=['.\save_image\',Date,'6 Mosaic of results','.jpg']; imwrite(I4,str);
    str='The results are saved in the save_image folder.\n\n'; fprintf(str);