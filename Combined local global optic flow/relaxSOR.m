
function [u,v,J,convergenceError] = relaxSOR(u, v, J, alpha, wFactor)
% 
%  relaxSOR
% 
%  SOR relaxation iteration for CLG-OF equations.
%  Each call to this function updates the current value of the solution,
%  u[1..m][1..n], v[1..m][1..n], using the motion tensor
%  J[JROWS][JCOLS][1..m][1..n].
%  Neumann boundary conditions are used (derivatives are set to zero).
% 
%  Parameters:
% 
%  u               Pointer to the optical flow horizontal component
%                  matrix/array.
%  v               Pointer to the optical flow vertical component matrix/array.
%  J[JROWS][JCOLS] Pointer to the array that stores the computed (and possibly
%                  smoothed) derivatives for each image/optical flow pixel.
%  nRows           Number of rows of the optical flow arrrays.
%  nCols           Number of columns of the optical flow arays.
%  alpha           Optical flow global smoothing coefficient.
%  wFactor         Relaxation parameter (if w=1.0 then the function is the
%                  Gauss-Seidel method).
% 

[nRows,nCols] = size(u);
error=0.0;

for i=2:nRows-1
    for j=2:nCols-1
        [u(i,j),v(i,j),errAt] = SOR_at(u, v, J, i, j, alpha, wFactor);
        error = error + errAt;
    end
end

[u,v,errBound] = boundaryCondition(u, v); 
error = error + errBound;

convergenceError = sqrt(error / (nRows*nCols));
end

function [ulocal,vlocal, error] = SOR_at(u,v,J,i,j, alpha, wFactor)


    %// h, could be less than 1.0
    h2a = 1.0/alpha;

    %// SOR formula
    numerator   = u(i,j-1) + u(i-1,j) + u(i,j+1) + u(i+1,j) - ...
                  h2a * (J(i,j,1,2) * v(i,j)+J(i,j,1,3));

    denominator = 4.0 + h2a * J(i,j,1,1);

    ulocal      = (1.0-wFactor) * u(i,j) + wFactor * numerator / denominator;

    numerator   = v(i,j-1) + v(i-1,j) + v(i,j+1) + v(i+1,j) - ...
                  h2a * (J(i,j,1,2) * ulocal+J(i,j,2,3));

    denominator = 4.0 + h2a * J(i,j,2,2);

    vlocal      = (1.0-wFactor) * v(i,j) + wFactor * numerator / denominator;

    error       = (ulocal - u(i,j)) * (ulocal - u(i,j));
    error      = error + (vlocal - v(i,j)) * (vlocal - v(i,j));

    u(i,j)     = ulocal;
    v(i,j)     = vlocal;

end

function [u,v,error] = boundaryCondition(u, v)

    [nRows,nCols] = size(u);
    error=0.0;

    % first and last rows
    for j=2:nCols-1

        i = 1;
        error = error + (u(i+1,j)-u(i,j)) * (u(i+1,j)-u(i,j));
        error = error + (v(i+1,j)-v(i,j)) * (v(i+1,j)-v(i,j));

        u(i,j) = u(i+1,j);
        v(i,j) = v(i+1,j);

        i = nRows;
        error = error + (u(i-1,j)-u(i,j)) * (u(i-1,j)-u(i,j));
        error = error + (v(i-1,j)-v(i,j)) * (v(i-1,j)-v(i,j));

        u(i,j) = u(i-1,j);
        v(i,j) = v(i-1,j);
    end

    % first and last columns
    for i=2:nRows-1

        j = 1;
        error = error + (u(i,j+1)-u(i,j)) * (u(i,j+1)-u(i,j));
        error = error + (v(i,j+1)-v(i,j)) * (v(i,j+1)-v(i,j));

        u(i,j) = u(i,j+1);
        v(i,j) = v(i,j+1);

        j = nCols;
        error = error + (u(i,j-1)-u(i,j)) * (u(i,j-1)-u(i,j));
        error = error + (v(i,j-1)-v(i,j)) * (v(i,j-1)-v(i,j));

        u(i,j) = u(i,j-1);
        v(i,j) = v(i,j-1);
    end


    % corners
    i=1; j=1;
    error = error + (u(i+1,j+1)-u(i,j)) * (u(i+1,j+1)-u(i,j));
    error = error + (v(i+1,j+1)-v(i,j)) * (v(i+1,j+1)-v(i,j));

    u(i,j) = u(i+1,j+1);
    v(i,j) = v(i+1,j+1);

    j = nCols;
    error = error + (u(i+1,j-1)-u(i,j)) * (u(i+1,j-1)-u(i,j));
    error = error + (v(i+1,j-1)-v(i,j)) * (v(i+1,j-1)-v(i,j));

    u(i,j) = u(i+1,j-1);
    v(i,j) = v(i+1,j-1);
    
    i = nRows;
    error = error + (u(i-1,j-1)-u(i,j)) * (u(i-1,j-1)-u(i,j));
    error = error + (v(i-1,j-1)-v(i,j)) * (v(i-1,j-1)-v(i,j));

    u(i,j) = u(i-1,j-1);
    v(i,j) = v(i-1,j-1);
    
    j=1;
    error = error + (u(i-1,j+1)-u(i,j)) * (u(i-1,j+1)-u(i,j));
    error = error + (v(i-1,j+1)-v(i,j)) * (v(i-1,j+1)-v(i,j));

    u(i,j) = u(i-1,j+1);
    v(i,j) = v(i-1,j+1);

end

