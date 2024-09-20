%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gao Chenzhong
% Contact: gao-pingqi@qq.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [magnitude,orientation] = W_Parital_Main_Orientation(I,R1,R2,s,int_flag)
sigma = 0.5;
w = 2*round(3*sigma)+1;
w = fspecial('gaussian',[w,w],sigma);
I = imfilter(I,w,'replicate');

Ns = 4; No = 6;
EO = LogGabor(I,Ns,No,3,1.6,0.75);  % minWaveLength = 3, mult = 1.6, sigmaOnf = 0.75
[M,N] = size(I); clear I
Gx = zeros(M,N); Gy = zeros(M,N);
angle = pi*(0:No-1)/No;
angle_cos = cos(angle);
angle_sin = sin(angle);
for j=1:No
    for i=1:Ns
        Gx = Gx - imag(EO{i,j}) * angle_cos(j);  % 注意图像与matlab矩阵坐标关系，y轴是反的
        Gy = Gy + imag(EO{i,j}) * angle_sin(j);
    end
end

W = floor(R2);  % 窗半径
dx = -W : W;  % 邻域x坐标
dy = -W : W;  % 邻域y坐标
[dx,dy] = meshgrid(dx,dy);
Wcircle = ((dx.^2 + dy.^2) < (W+1)^2)*1.0;  % 圆形窗
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

Gxx = imfilter(Gx.*Gx, h, 'replicate');
Gyy = imfilter(Gy.*Gy, h, 'replicate');
Gxy = imfilter(Gx.*Gy, h, 'replicate'); clear Gx Gy
Gsx = Gxx-Gyy;                          clear Gxx Gyy
Gsy = 2*Gxy;                            clear Gxy

orientation = atan2(Gsy,Gsx)/2 + pi/2;  % 取值范围：[-pi,pi] ——> [0,pi]
orientation = mod(orientation,pi);  % 取值范围：[0,pi] ——> [0,pi)
if int_flag
    magnitude = ones(M,N);
else
    magnitude = sqrt(sqrt(Gsx.^2+Gsy.^2));
end