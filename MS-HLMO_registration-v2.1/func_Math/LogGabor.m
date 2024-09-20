%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gao Chenzhong
% Contact: gao-pingqi@qq.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EO = LogGabor(img, nscale, norient, minWaveLength, mult, sigmaf)
% default parameters:
% nscale        = 4;     % Number of wavelet scales.
% norient       = 6;     % Number of filter orientations.
% minWaveLength = 3;     % Wavelength of smallest scale filter.
% mult          = 2.1;   % Scaling factor between successive filters.
% sigmaf        = 0.55;  % Ratio of the standard deviation of the
                         % Gaussian describing the logGabor filter's
                         % transfer function in the frequency domain
                         % to the filter center frequency.

[rows,cols,~] = size(img);
imagefft = fft2(img); clear img
EO = cell(nscale, norient);

if mod(cols,2)
    xrange = (-(cols-1)/2:(cols-1)/2)/(cols-1);
else
    xrange = (-cols/2:(cols/2-1))/cols;
end

if mod(rows,2)
    yrange = (-(rows-1)/2:(rows-1)/2)/(rows-1);
else
    yrange = (-rows/2:(rows/2-1))/rows;
end

[x,y] = meshgrid(xrange, yrange); clear xrange yrange
radius = sqrt(x.^2 + y.^2);
theta = atan2(-y,x); clear x y

radius = ifftshift(radius);
theta  = ifftshift(theta);
radius(1,1) = 1;  % log friendly

sintheta = sin(theta);
costheta = cos(theta);
clear theta

lp = lowpassfilter([rows,cols],0.45,15);  % Radius 0.45, 'sharpness' 15

logGabor = cell(1,nscale);
for s = 1:nscale
    wavelength = minWaveLength*mult^(s-1);
    fo = 1.0/wavelength;
    logGabor{s} = exp((-(log(radius/fo)).^2) / (2 * log(sigmaf)^2));
    logGabor{s} = logGabor{s}.*lp;
    logGabor{s}(1,1) = 0;  % log friendly
end
clear radius lp

for o = 1:norient
    angle = (o-1)*pi/norient;
    ds = sintheta * cos(angle) - costheta * sin(angle);
    dc = costheta * cos(angle) + sintheta * sin(angle);
    dtheta = abs(atan2(ds,dc));
    dtheta = min(dtheta*norient/2,pi);
    spread = (cos(dtheta)+1)/2;
    for s = 1:nscale
        filter = logGabor{s} .* spread;
        EO{s,o} = ifft2(imagefft .* filter);
    end
end