%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gao Chenzhong
% Contact: gao-pingqi@qq.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EO = LogGabor(varargin)
    
% Get arguments and/or default values
[im, nscale, norient, minWaveLength, mult, sigmaOnf] = checkargs(varargin(:));

[rows,cols] = size(im);
imagefft = fft2(im);              % Fourier transform of image
EO = cell(nscale, norient);       % Array of convolution results
% Pre-compute some stuff to speed up filter construction

% Set up X and Y matrices with ranges normalised to +/- 0.5
% The following code adjusts things appropriately for odd and even values
% of rows and columns.
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

[x,y] = meshgrid(xrange, yrange);

radius = sqrt(x.^2 + y.^2);       % Matrix values contain *normalised* radius from centre.
theta = atan2(-y,x);              % Matrix values contain polar angle.
                                  % (note -ve y is used to give +ve
                                  % anti-clockwise angles)

radius = ifftshift(radius);       % Quadrant shift radius and theta so that filters
theta  = ifftshift(theta);        % are constructed with 0 frequency at the corners.
radius(1,1) = 1;                  % Get rid of the 0 radius value at the 0
                                  % frequency point (now at top-left corner)
                                  % so that taking the log of the radius will
                                  % not cause trouble.
sintheta = sin(theta);
costheta = cos(theta);
clear x; clear y; clear theta;    % save a little memory

% Filters are constructed in terms of two components.
% 1) The radial component, which controls the frequency band that the filter
%    responds to
% 2) The angular component, which controls the orientation that the filter
%    responds to.
% The two components are multiplied together to construct the overall filter.

% Construct the radial filter components...
% First construct a low-pass filter that is as large as possible, yet falls
% away to zero at the boundaries.  All log Gabor filters are multiplied by
% this to ensure no extra frequencies at the 'corners' of the FFT are
% incorporated as this seems to upset the normalisation process when
% calculating phase congrunecy.
lp = lowpassfilter([rows,cols],.45,15);   % Radius .45, 'sharpness' 15

logGabor = cell(1,nscale);

for s = 1:nscale
    wavelength = minWaveLength*mult^(s-1);
    fo = 1.0/wavelength;                  % Centre frequency of filter.
    logGabor{s} = exp((-(log(radius/fo)).^2) / (2 * log(sigmaOnf)^2));
    logGabor{s} = logGabor{s}.*lp;        % Apply low-pass filter
    logGabor{s}(1,1) = 0;                 % Set the value at the 0 frequency point of the filter
                                          % back to zero (undo the radius fudge).
end

%% The main loop...
for o = 1:norient                    % For each orientation...
    % Construct the angular filter spread function
    angl = (o-1)*pi/norient;           % Filter angle.
    % For each point in the filter matrix calculate the angular distance from
    % the specified filter orientation.  To overcome the angular wrap-around
    % problem sine difference and cosine difference values are first computed
    % and then the atan2 function is used to determine angular distance.
    ds = sintheta * cos(angl) - costheta * sin(angl);    % Difference in sine.
    dc = costheta * cos(angl) + sintheta * sin(angl);    % Difference in cosine.
    dtheta = abs(atan2(ds,dc));                          % Absolute angular distance.
    % Scale theta so that cosine spread function has the right wavelength and clamp to pi
    dtheta = min(dtheta*norient/2,pi);
    % The spread function is cos(dtheta) between -pi and pi.  We add 1,
    % and then divide by 2 so that the value ranges 0-1
    spread = (cos(dtheta)+1)/2;

    for s = 1:nscale                  % For each scale...
        filter = logGabor{s} .* spread;      % Multiply radial and angular
                                             % components to get the filter.
        % Convolve image with even and odd filters returning the result in EO
        EO{s,o} = ifft2(imagefft .* filter);
    end
end  % For each orientation


%%------------------------------------------------------------------
% CHECKARGS
%
% Function to process the arguments that have been supplied,
% assign default values as needed and perform basic checks.
    
function [im, nscale, norient, minWaveLength, mult, sigmaOnf] = checkargs(arg)

nargs = length(arg);

if nargs < 1
    error('No image supplied as an argument');
end    

% Set up default values for all arguments and then overwrite them
% with with any new values that may be supplied
im              = [];
nscale          = 4;     % Number of wavelet scales.
norient         = 6;     % Number of filter orientations.
minWaveLength   = 3;     % Wavelength of smallest scale filter.
mult            = 2.1;   % Scaling factor between successive filters.
sigmaOnf        = 0.55;  % Ratio of the standard deviation of the
                         % Gaussian describing the log Gabor filter's
                         % transfer function in the frequency domain
                         % to the filter center frequency.

% Allowed argument reading states
allnumeric   = 1;       % Numeric argument values in predefined order
keywordvalue = 2;       % Arguments in the form of string keyword
                        % followed by numeric value
readstate = allnumeric; % Start in the allnumeric state

if readstate == allnumeric
    for n = 1:nargs
        if isa(arg{n},'char')
            readstate = keywordvalue;
            break;
        else
            if     n == 1, im            = arg{n};
            elseif n == 2, nscale        = arg{n};
            elseif n == 3, norient       = arg{n};
            elseif n == 4, minWaveLength = arg{n};
            elseif n == 5, mult          = arg{n};
            elseif n == 6, sigmaOnf      = arg{n};
            end
        end
    end
end

% Code to handle parameter name - value pairs
if readstate == keywordvalue
    while n < nargs

        if ~isa(arg{n},'char') || ~isa(arg{n+1}, 'double')
            error('There should be a parameter name - value pair');
        end

        if     strncmpi(arg{n},'im'      ,2), im =        arg{n+1};
        elseif strncmpi(arg{n},'nscale'  ,2), nscale =    arg{n+1};
        elseif strncmpi(arg{n},'norient' ,4), norient =   arg{n+1};
        elseif strncmpi(arg{n},'minWaveLength',2), minWaveLength = arg{n+1};
        elseif strncmpi(arg{n},'mult'    ,2), mult =      arg{n+1};
        elseif strncmpi(arg{n},'sigmaOnf',2), sigmaOnf =  arg{n+1};
        else, error('Unrecognised parameter name');
        end

        n = n+2;
        if n == nargs
            error('Unmatched parameter name - value pair');
        end

    end
end

if isempty(im)
    error('No image argument supplied');
end

if ndims(im) == 3
    warning('Colour image supplied: converting image to greyscale...')
    im = double(rgb2gray(im));
end

if ~isa(im, 'double')
    im = double(im);
end

if nscale < 1
    error('nscale must be an integer >= 1');
end

if norient < 1 
    error('norient must be an integer >= 1');
end

if minWaveLength < 2
    error('It makes little sense to have a wavelength < 2');
end