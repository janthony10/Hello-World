
function [res,Uout,Vout] = calcCLG_of(prevFrame,currFrame,para)
    
DEFAULT_GAUSSIAN_WINDOW_SIZE = 2.5;
MIN_FILTER_SIZE_RHO = 3;
MIN_FILTER_SIZE_SIGMA = 3;
EPS = 1E-12;
MIN_ERROR = 0.0001;

rho = para.rho;
iterations = para.numIt;
alpha = para.alpha;
wFactor = para.w;
verbose = para.verbose;
coupledMode = para.coupledMode;

    [nRows,nCols] = size(prevFrame);
    if (verbose)
        fprintf('calc_clg\n');
        fprintf('  setting up variables\n');
    end

    % h, could be less than 1.0;
    h=1.0;

    u = zeros(nRows,nCols);
    v = zeros(nRows,nCols);
    Uout = u;
    Vout = v;
    
       
    % J coeffs obtention
    J = computeDerivatives(prevFrame, currFrame, verbose);

    if (verbose)
        fprintf('local spatio temporal smoothing\n');
    end

    if (rho > 0)
        % odd filter size
        filterSize = 2*uint8(DEFAULT_GAUSSIAN_WINDOW_SIZE * rho) + 1;
        
        if (verbose)
            fprintf('  filter size %i for rho %f\n', filterSize, rho);
        end

        if (filterSize < MIN_FILTER_SIZE_RHO)
            filterSize = MIN_FILTER_SIZE_RHO;
            if (verbose)
                fprintf('  after checking, new filter size %i for rho %f\n',filterSize, rho);
            end
        end
                
        for k = 1:9
            J(:,:,k) = gaussian2(J(:,:,k),rho);
        end
        
    end

    if (iterations == 0)
        iterations = uint8(nRows * nCols / 8.0);
    end

    if (verbose)
        fprintf('  performing %i relax iterations\n', iterations);
    end

    count = 0;
    error = 1.0000;
    convergenceError = 0.0;
    
    while (count<iterations && convergenceError*1.01<error)

        if (count > 0)
            error = convergenceError;
        end
        
        [u,v,J,convergenceError] = relaxSOR(u, v, J, alpha, wFactor);
        
        count = count + 1; 
    end

    %% Fill output values.
    if (verbose)
        fprintf('  filling output after %d iterations, error=%f\n', count, error);
    end
    %% the sign change is because the y-axis is inverted
    for i=1:nRows
        for j=1:nCols
            Uout(i,j) = Uout(i,j) + u(i,j);
            Vout(i,j) = Vout(i,j) + v(i,j);
        end
    end
   
    if (verbose)
        fprintf('calc_clg: done\n');
    end
    res = 1;

    end

    function J = computeDerivatives(prevFrame,currFrame,verbose)
    [nRows,nCols] = size(prevFrame);
    J = zeros(nRows, nCols, 3, 3);
    
    h = 1;
    k = [-1 8 0 -8 1]/12*h;
    
    Fx = conv2(currFrame,k,'same');
    Fy = conv2(currFrame,k,'same');
    Ft = currFrame - prevFrame;
    
    if (verbose)
        fprintf('    storing J coeffs.\n');
    end
    
    J(:,:,1,1) = Fx .* Fx;
    J(:,:,1,2) = Fx .* Fy;
    J(:,:,1,3) = Fx .* Ft;
    J(:,:,2,1) = J(:,:,1,2);
    J(:,:,2,2) = Fy .* Fy;
    J(:,:,2,3) = Fy .* Ft;
    J(:,:,3,1) = J(:,:,1,3);
    J(:,:,3,2) = J(:,:,2,3);
    J(:,:,3,3) = Ft .* Ft;
    
    if (verbose)
        fprintf('    done\n');
    end
    
    end
   
