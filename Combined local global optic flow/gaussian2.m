function img_smooth = gaussian2(I1, sigma)
DEFAULT_GAUSSIAN_WINDOW_SIZE = 2.5;

if size(I1,3) ~= 1
    error('The input image should be in grayscale!');
end

window_size = DEFAULT_GAUSSIAN_WINDOW_SIZE;
Fsize = uint8(window_size * sigma) + 1;
den  = 2*sigma*sigma;

x = double(-Fsize:1:Fsize);

B = 1 / (sigma * sqrt(2.0 * pi)) * exp(-x .* x / den);

B = B / sum(B);

img_smooth=conv2(B,B,I1,'same');





end