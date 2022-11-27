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

%% Read Images
[image_1, image_2] = Readimage;

%% Image Preproscessing
resample1 = 1/1; resample2 = 1/1;
[I1_s,I1,I2_s,I2] = Preproscessing(image_1,resample1,image_2,resample2);
figure,imshow(I1_s),title('Reference image'); pause(0.01)
figure,imshow(I2_s),title('Sensed Image'); pause(0.01)
% figure,imshow(I1),title('Reference image base'); pause(0.01)
% figure,imshow(I2),title('Sensed Image base'); pause(0.01)

%%
warning off
    fprintf('\n** Registration starts, have fun\n\n'); ts=cputime;

[cor1,cor2] = main_image_matching(I1_s,I1,I2_s,I2,...
    int_flag,rot_flag,scl_flag,par_flag,key_type);
% matchment0 = Show_Matches(I1,I2,cor1,cor2,1);pause(0.01)
cor1 = cor1/resample1; cor2 = cor2/resample2;
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
        str=['.\save_image\',Date,'3 Reference Image','.jpg']; imwrite(I1_s,str);
    otherwise
        str=['.\save_image\',Date,'1 Reference Image','.mat']; save(str,'I1_r');
        str=['.\save_image\',Date,'3 Reference Image','.jpg']; imwrite(I1_rs,str);
end
str=['.\save_image\',Date,'2 Registered Image','.mat']; save(str,'I2_r');
str=['.\save_image\',Date,'4 Registered Image','.jpg']; imwrite(I2_rs,str);
str=['.\save_image\',Date,'5 Fusion of results','.jpg']; imwrite(I3,str);
str=['.\save_image\',Date,'6 Mosaic of results','.jpg']; imwrite(I4,str);
    str='The results are saved in the save_image folder.\n\n'; fprintf(str);