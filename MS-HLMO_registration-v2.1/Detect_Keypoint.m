%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gao Chenzhong
% Contact: gao-pingqi@qq.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kps = Detect_Keypoint(I,scale,thresh,radius,N,nOctaves,G_resize,type,disp)
if disp
    I_o = I;
end

%% Phase-Congruency participate
if contains(type, 'PC')
    [pc_M,pc_m,~,~,~,~,~] = phasecong3(I, 'nscale',4, 'norient',6,...
        'minWaveLength',3, 'mult',1.6, 'sigmaOnf',0.75, 'g', 3, 'k',1);
    I = pc_M + pc_m;
end

%% Image norm and mask
I = I - min(I(:));
I = (I/mean(I(:))/2)*255;
msk = Mask(I,-10);

%% Keypoints detectation
if contains(lower(type), 'harris')
    kps = Harris(I,scale,thresh,radius);
elseif contains(lower(type), 'shitomasi')
    kps = ShiTomasi(I,scale,thresh,radius);
else
    assert(false,'Unexpected Feature Point Type encountered.');
end

%% Post-processing
kps = Remove_Boundary_Points(kps,msk,max(scale,G_resize^(nOctaves-2)));
if size(kps,1)<10
    kps = []; return
end
kps = sortrows(kps,-3);
kps = kps(1:min(N,size(kps,1)),:);

%% Show detected keypoints
if disp==1
    figure, imshow(I_o,[]), hold on, plot(kps(:,1),kps(:,2),'r.');
%     for i=1:size(kps,1)
%         text(kps(i,1),kps(i,2),num2str(i),'color','y');
%     end
    title(['Detected keypoints: ',num2str(size(kps,1))]); drawnow
end


function msk = Mask(I,th)
I = I./max(I(:))*255;
msk = double(I>th);
h = D2gauss(7,4,7,4,0);
msk = (conv2(msk,h,'same')>0.0);  % default:0.8


function p = Remove_Boundary_Points(loc,msk,s)
se = strel('disk',s);
msk = ~(imdilate(~msk,se));
p = [];
for i = 1:size(loc,1)
    if msk(loc(i,2),loc(i,1)) == 1
        p = [p;loc(i,:)];
    end
end


function I_p = Image_Pat(I,s)
[m,n] = size(I);
I_p = zeros([m+2*s,n+2*s]);
I_p(s+1:end-s,s+1:end-s) = I;