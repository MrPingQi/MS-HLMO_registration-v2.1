%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gao Chenzhong
% Contact: gao-pingqi@qq.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gao Chenzhong
% Contact: gao-pingqi@qq.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function descriptors = Multiscale_Descriptor(I,kps,patch_size,NBA,NBO,...
    nOctaves,nLayers,G_resize,G_sigma,int_flag,rot_flag,par_flag)

sig = Get_Gaussian_Scale(G_sigma,nLayers);

descriptors = cell(nOctaves,nLayers);
if par_flag
    nScales = nOctaves*nLayers;
    descriptor = cell(nScales,1);
    parfor k=1:nScales
        octave = ceil(k/nLayers);
        layer = k-(octave-1)*nLayers;
        if octave>1
            kps_t = round(kps(:,1:2)./G_resize^(octave-1));
            [~,index,~] = unique(kps_t,'rows');
        else
            kps_t = kps(:,1:2);
            index = 1:size(kps,1);
        end
        I_t = Gaussian_Scaling_par(I, octave, layer, G_resize, sig(layer));
        descriptor{k} = GGLOHv2_Descriptor(I_t, ...
            [kps_t(index,1:2), kps(index,1:2)],patch_size,NBA,NBO,int_flag,rot_flag);
    end
    for k=1:nScales
        octave = ceil(k/nLayers);
        layer = k-(octave-1)*nLayers;
        descriptors{octave,layer} = descriptor{k};
    end

else
    I_t = [];
    for octave=1:nOctaves
        kps_t = round(kps(:,1:2)./G_resize^(octave-1));
        [~,index,~] = unique(kps_t,'rows');
        for layer=1:nLayers
            I_t = Gaussian_Scaling(I, I_t, octave, layer, G_resize, sig(layer));
            descriptors{octave,layer} = GGLOHv2_Descriptor(I_t, ...
                [kps_t(index,1:2), kps(index,1:2)],patch_size,NBA,NBO,int_flag,rot_flag);
        end
    end
end


function sig = Get_Gaussian_Scale(sigma,numLayers)
sig = zeros(1,numLayers);
sig(1) = sigma;  % 认为第一个图像尺度就是σ
if numLayers<2
    return
end
k = 2^(1.0/(numLayers-1));
for i = 2:1:numLayers
    sig_prev = k^(i-2)*sigma;
    sig_curr = k*sig_prev;
    sig(i) = sqrt(sig_curr^2-sig_prev^2);
end


function I_t = Gaussian_Scaling(I,I_t,Octave,Layer,G_resize,sig)
if(Octave==1 && Layer==1)
    I_t = I;
elseif(Layer==1)
    I_t = imresize(I,1/G_resize^(Octave-1),'bicubic');
else
    window_gaussian = round(2*sig);
    window_gaussian = 2*window_gaussian+1;
    w = fspecial('gaussian',[window_gaussian,window_gaussian],sig);
    I_t = imfilter(I_t,w,'replicate');
end


function I = Gaussian_Scaling_par(I,Octave,Layer,G_resize,sig)
if(Octave==1 && Layer==1)
    return
elseif(Octave==1)
    window_gaussian = round(2*sig);
    window_gaussian = 2*window_gaussian+1;
    w = fspecial('gaussian',[window_gaussian,window_gaussian],sig);
    I = imfilter(I,w,'replicate');
elseif(Layer==1)
    I = imresize(I,1/G_resize^(Octave-1),'bicubic');
else
    I = imresize(I,1/G_resize^(Octave-1),'bicubic');
    window_gaussian = round(2*sig);
    window_gaussian = 2*window_gaussian+1;
    w = fspecial('gaussian',[window_gaussian,window_gaussian],sig);
    I = imfilter(I,w,'replicate');
end