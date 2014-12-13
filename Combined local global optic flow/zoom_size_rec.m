function [nx,ny] = zoom_size_rec(nx,ny,factor,iter)

for i = 1:iter
    [nxx,nyy] = zoom_size(nx,ny,factor);
    nx = nxx;
    ny = nyy;
end


end