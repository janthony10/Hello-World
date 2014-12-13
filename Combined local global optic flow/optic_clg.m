function [u,v] = optic_clg(im1,im2,para)

if nargin > 3
    error('There should be three parameters.')
end

if nargin == 2
    para.alpha = 5;             %regularization parameter
    para.rho = 1;               %local smoothing parameter
    para.sigma = 1.4;           %pre-processing smoothing parameter
    para.numIt = 1000;          %Iterations numbers
    para.w = 1.9;               %SOR relaxation factor, between 0 and 2. Default: 1.9.
    para.nScales = 1;           %number of scales for the pyramidal approach. Default: 1.
    para.scaleFactor = 0.5;     %scale factor, between 0 and 1. Default: 0.5.
    para.coupledMode = 0;       %iteration type: 1->Gauss-Seidel, 0->SOR. Default: 1.
    para.verbose = 1;           %shows (1) or hide (0) messages. Default: 1.
end

im1 = im2double(im1);
im2 = im2double(im2);

[h,w,~] = size(im1);

% at the highest pyramid level the image must be have
% at least 45 pixels side

max_n_scales = uint8(floor(1.0 + log(30.0 / min(1.0*w, 1.0*h)) / log(para.scaleFactor)));

if (para.nScales > max_n_scales)
    para.nScales = max_n_scales;
    if (verbose)
        printf('nScales corrected to the max value allowed: %d\n', para.nScales);
    end
end

[res,u,v] = calcMSCLG_OF(im1, im2, para);

end