global mkl gpu;

startup_STB();
% Parameters------------------------------
Fr = 0.9;               % The Froude number
source = [0,0,-1];      % The location of the source [x,y,z]
epsilon = [0.02:0.02:0.92];  % Epsilon can be a column vector where the solution to epsilon(i-1) is used as the initial guess for epsilon(i)
singType = false;           % The type of singularity, true = source, false = dipole
intMethod = [];         % The weighting scheme for numerical integration
intMethod.x = 'trap';
intMethod.y = 'trap';
n = 0.005;              % The upstraem dampening parameter

% Toggle these variables depending if you have Intel MKL or a gpu that
% can run CUDA 5.0
mkl = true;
gpu = true;


%% Create the mesh-----------------------------
% Intel MKL is required to use a banded storage preconditioner, otherwise
% you must use a dense storage preconditioner which will restrict the size
% of the mesh that can be used.
if mkl
    N = 735;
    M = 237;
    x0=-18;
    xEnd = 32.4;
    deltaX = (xEnd-x0)/(N-1);
    deltaY = deltaX;
    band = 39;
    preconType = 'band';
    
    % This code selects the largest bandwidth posible with the given amount
    % of system memory (only possible with a 'band' preconditioner)
%     preconditionerMaxMem = 85;
%     bandf = @(b) calcMem(M,N,b,'band')-preconditionerMaxMem;
%     band = floor(fzero(bandf,[0,M]));
else
    error('Intel MKL must be used');
    
end

% Display the memory used by the preconditioner
mem = calcMem(M,N,band,preconType); 
disp(['Approx mem used by preconditioner ',num2str(mem),'Gbs']);

% Flat surface initial guess
yDash0 = [repmat([x0;ones(N,1)],M,1);zeros(M*(N+1),1)];

% The folder that the surface data is saved to
folder = './Stats/Dipole-Fr09/';
% Calculate solutions
computeSurface(M,N,deltaX,deltaY,x0,source,epsilon,singType,Fr,intMethod,n,band,preconType,yDash0,folder);

