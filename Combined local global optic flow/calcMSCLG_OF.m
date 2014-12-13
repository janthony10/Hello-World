function [status,us_hr,vs_hr] = calcMSCLG_OF(image1, image2, para)

status = 0;
[nRows,nCols,~] = size(image1);

DEFAULT_GAUSSIAN_WINDOW_SIZE = 2.5;
MIN_FILTER_SIZE_RHO = 3;
MIN_FILTER_SIZE_SIGMA = 3;
EPS = 1E-12;
MIN_ERROR = 0.0001;
iterations = para.numIt;
alpha = para.alpha;
rho = para.rho;
sigma = para.sigma;
wFactor = para.w;
nScales = para.nScales;
scaleFactor = para.scaleFactor;
coupledMode = para.coupledMode;
verbose = para.verbose;

if (verbose)
    fprintf('calcMSCLG_OF\n');
    fprintf('  parameters: nScales=%d, scaleFactor=%f\n', nScales, scaleFactor);
end

%% variables declaration

Isize = nCols * nRows;
I1s = image1;
I1s_hr = I1s;
I2s = image2;
I2s_hr = I2s;

% us = uOut;
% us_hr = us;
% vs = vOut;
% us_hr = vs;
nxx = zeros(nScales,1);
nyy = zeros(nScales,1);

nxx(1) = nCols;
nyy(1) = nRows;



%% pre-smoothing the finest scale images
filterSize = uint8(2*(DEFAULT_GAUSSIAN_WINDOW_SIZE * sigma) + 1);

if (filterSize >= MIN_FILTER_SIZE_SIGMA)
    fprintf('  apply gaussian smoothing for each frame\n');
    I1s = gaussian2(I1s,sigma);
    I2s = gaussian2(I2s,sigma);
    
else
    fprintf('  no gaussian smoothing was applied to input frames\n');
end

%% create the scales
for (s=2:nScales)
    
    if (verbose)
        fprintf('  scales %d init\n', s);
    end
    
    [nxx(s),nyy(s)] = zoom_size(nxx(s-1), nyy(s-1), scaleFactor);
    
    
    % compute the zoom from the previous finer scale
%     I1s_hr = zoom_out(I1s, scaleFactor);
%     I2s_hr = zoom_out(I2s, scaleFactor);
end
% initialize the flow

%     for (i=0; i < nxx[nScales-1] * nyy[nScales-1]; i++) {
%         us[nScales-1][i] = 0.0;
%         vs[nScales-1][i] = 0.0;
%     }
%
%% pyramidal approximation to the optic flow

for s=nScales:-1:1
%     [nxx_hr,nyy_hr] = zoom_size_rec(nxx, nyy, scaleFactor,s);
    if (verbose)
        fprintf('Scale: %d %dx%d, nScales=%d\n', s, nxx(s), nyy(s), nScales);
    end
    
    I1s_hr = zoom_out(I1s,[nyy(s) nxx(s)]);
    I2s_hr = zoom_out(I2s,[nyy(s) nxx(s)]);
    if s == nScales
        us_hr = zeros(nyy(s),nxx(s));
        vs_hr = zeros(nyy(s),nxx(s));
    end
    % warp the second image before computing derivatives
    image2warped = imWarp(I2s_hr,us_hr,vs_hr);
    
    % compute the optical flow
    [~,us_hr,vs_hr] = calcCLG_of(I1s_hr,image2warped,para);
    
    % if this was the last scale, finish now
    if s == 1;
        break;
    end
    
    % otherwise, upsample the optical flow
    % zoom the optic flow for the next finer scale
    
    us = zoom_in(us_hr,[nyy(s-1) nxx(s-1)]);
    vs = zoom_in(vs_hr,[nyy(s-1) nxx(s-1)]);
    
    
    % scale the optic flow with the appropriate zoom factor
    us = us * 1.0 / scaleFactor;
    vs = vs * 1.0 / scaleFactor;
    
    
    if (verbose)
        fprintf('calcMSCLG_OF: done\n');
    end
    
    us_hr = us;
    vs_hr = vs;
end

status = 1;
end

function [nxx,nyy] = zoom_size_rec(nx,ny,factor,iter)

for i = 1:iter
    [nx,ny] = zoom_size(nx,ny,factor);
end

nxx = nx;
nyy = ny;

end

function[nxx,nyy] = zoom_size(nx, ny,factor)
%Compute the size of a zoomed image from the zoom factor

nxx = floor(nx * factor);
nyy = floor(ny * factor);

end

function Iout = zoom_out_rec(I, factor,iter)

for i = 1:iter
    I = zoom_out(I,factor);
end
Iout = I;

end

function Iout = zoom_out(I, new_size)
%Downsample an image
ZOOM_SIGMA_ZERO = 0.6;
[ny,nx] = size(I);
% compute the size of the zoomed image
% [nxx,nyy] = zoom_size(nx, ny, factor);

factor =mean([new_size(1)/ny, new_size(2)/nx]);
% compute the Gaussian sigma for smoothing
sigma = ZOOM_SIGMA_ZERO * sqrt(1.0/(factor*factor) - 1.0);

% pre-smooth the image
I = gaussian2(I,sigma);

% re-sample the image using bicubic interpolation
%#pragma omp parallel for

Iout = imresize(I,new_size,'bicubic');

end

function Iout = zoom_in(I, new_size)

% int nxx,        // width of the zoomed image
% int nyy         // height of the zoomed image

[ny,nx] = size(I);
% nxx = round(nx/factor);
% nyy = round(ny/factor);


% re-sample the image using bicubic interpolation
% #pragma omp parallel for

Iout = imresize(I,new_size, 'bicubic');



end