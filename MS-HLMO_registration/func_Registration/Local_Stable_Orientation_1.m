%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gao Chenzhong
% Contact: gao-pingqi@qq.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [magnitude,orientation] = Local_Stable_Orientation_1(I,R1,R2,s)

sigma=0.5;
w=2*round(3*sigma)+1;
w=fspecial('gaussian',[w,w],sigma);
I=imfilter(I,w,'replicate');

Ns=4; No=6;
eo = LogGabor(I,Ns,No,3,'mult',1.6,'sigmaOnf',0.75);
[M,N] = size(I);
Gx = zeros(M,N); Gy = zeros(M,N);
for j=1:No
    for i=1:Ns
        Gx = Gx-imag(eo{i,j})*cos(pi*(j-1)/No); % 注意图像与matlab矩阵坐标关系，y轴是反的
        Gy = Gy+imag(eo{i,j})*sin(pi*(j-1)/No);
    end
end

W = floor(R2); % 窗半径
dx = -W : W; % 邻域x坐标
dy = -W : W; % 邻域y坐标
[dx,dy] = meshgrid(dx,dy);
Wcircle = ((dx.^2 + dy.^2) < (W+1)^2)*1.0; % 圆形窗
Patchsize = 2*W+1;

if s==1
    h = fspecial('gaussian',[Patchsize,Patchsize], R1/6);
else
    step = (R2-R1)/(s-1);
    h = zeros(Patchsize,Patchsize);
    for i=0:s-1
        sigma = (R1+step*i)/6;
        h = h + fspecial('gaussian',[Patchsize,Patchsize], sigma);
    end
end
h = h.*Wcircle;

Gxx = conv2(Gx.*Gx, h, 'same');
Gyy = conv2(Gy.*Gy, h, 'same');
Gxy = conv2(Gx.*Gy, h, 'same');
Gsx = Gxx-Gyy; Gsy = 2*Gxy;

magnitude = sqrt(sqrt(Gsx.^2+Gsy.^2));
orientation = atan2(Gsy,Gsx)/2 + pi/2; % 取值范围：[-pi,pi] ——> [0,pi]