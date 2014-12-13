function[nxx,nyy] = zoom_size(nx, ny,factor)
%Compute the size of a zoomed image from the zoom factor

nxx = round(nx * factor + 0.5);
nyy = round(ny * factor + 0.5);

end