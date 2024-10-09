function yDash = computeSurface(M,N,deltaX,deltaY,x0,source,epsilon,singType,Fr,intMethod,n,band,preconType,yDash0,folder)
% computeSurface
% M,N - Dimensions of the mesh
% deltaX, deltaY - The x,y spacing of the mesh
% x0 - The location of the upstream point
% source - The location of the source [x,y,z]
% epsilon - The nonlinearity parameter epsilon
% singType - (bool) The type of singularity, true=source, false=dipole
% Fr - The Froude number
% intMethod - (struct) .x .y the weighting schemes for numerical
% integration in the x and y direction - values 'trap', 'simpson', 'boole'
% n - upstream dampening parameter
% band - the block-bandwidth of the preconditiontioner values [0,M-1]
% preconType - The type of preconditioner used - values 'band' (banded
% storage preconditioner, requires Intel MKL), 'dense' (dense storage) and
% 'full' (full Newton's method)
% yDash0 - The initial guess
% folder - The folder that the surface data is saved to

global PreSolveTime LinPreRunTime NonLinFuncRunCount NonLinFuncRunTime data nonLin kern baseMult mkl gpu base;
setTo0(true);
nonLin = 0;

% Input check
if band>(M-1)
    band=M-1;
end

if ~mkl
    if strcmp(preconType,'band')
        error('Intel MKL required for the banded storage preconditioner');
    end
end


% Initialise Function---------------------
if gpu
    kern = parallel.gpu.CUDAKernel('fGutsGPU.ptx','fGutsGPU.cu','guts');
    kern.ThreadBlockSize = [512,1,1];               
    kern.GridSize = [N-1,M,1];                                          
    kern.SharedMemorySize = 10000;

    kernbaseMult = parallel.gpu.CUDAKernel('baseMult.ptx','baseMult.cu','guts');
    kernbaseMult.ThreadBlockSize = [256,1,1];               
    kernbaseMult.GridSize = [N-1,M,1];                                   
    kernbaseMult.SharedMemorySize = 10000;

    baseMult = @(r) gather(feval(kernbaseMult,N*ones(1,1,'int64'),M*ones(1,1,'int64'),deltaX,r,zeros(size(r))));
else
    base = tril(ones(N-1)-3/4*eye(N-1)-1/4*diag(ones(N-2,1),-1));
    base = deltaX*[[1/4;1/2*ones(N-2,1)],base];
    base = [zeros(2,N+1);ones(N-1,1),base];
    baseMult = @(r) baseMultCPU(r);
end

psetfnBlock = @ (y,yscale,Fy,fscale) LinearBlockJacPre(y,Fy,M,N,deltaX,deltaY,x0,Fr,intMethod,n,yscale,fscale,band);
psetfnBand = @ (y,yscale,Fy,fscale) LinearBandJacPre(y,Fy,M,N,deltaX,deltaY,x0,Fr,intMethod,n,yscale,fscale,band);


%%
data = [];
data.kl = [];
data.A = [];
data.ipvtA = [];
data.Bcoef = [];
data.Ccoef = [];
data.D = [];
data.ipvtD = [];

%number of equations
neq = 2*(N+1)*M;

% Maximum Krylov subspace size
maxl = 100;

% Maximum norm of the search direction
mxnewt = 5;

% Scaling vectors
yscale = ones(neq,1);
yscale(1:(N+1):(M*(N+1))) = 1/abs(x0);
fscale = ones(neq,1);

strategy = 'LineSearch';

%%

% Set global values to zero
setTo0(true);

% Count the number of bootstrapping steps
[~,numEps] =size(epsilon);

% Generate the preconditioner
if strcmp(preconType,'dense')
    psetfnBlock(0,0,0,0);
elseif strcmp(preconType,'band')
    psetfnBand(0,0,0,0);
end

% Create folder to store solutions
s=mkdir(folder);
% itterate over all bootstrapping steps
yDash=yDash0;
for i=1:numEps
    
    % Set global values (not includin preconditioner set up time) to zero
    setTo0(false);
    JFNKLinBandTime = tic;
    eps = epsilon(:,i);
    f = @(x) BIFunction(x,source,eps,singType,M,N,deltaX,deltaY,x0,Fr,intMethod,n);
    
    % Set up KINSol options depeding on the preconditioner storage type
    if strcmp(preconType,'dense')   % Dense storage
        optionsFull = KINSetOptions('Verbose',true, ...
            'MaxNewtonStep',mxnewt, ...
                                'LinearSolver', 'GMRES', ...
                                'KrylovMaxDim', maxl,...
                                'MaxNumSetups',1,...
                                'PrecSolveFn',@preSolveBlock);
    elseif strcmp(preconType,'band')    %Banded storage
        optionsFull = KINSetOptions('Verbose',true, ...
            'MaxNewtonStep',mxnewt, ...
                                'LinearSolver', 'GMRES', ...
                                'KrylovMaxDim', maxl,...
                                'MaxNumSetups',1,...
                                'PrecSolveFn',@preSolveBnd);
    else                                % Finite difference Jacobian
        optionsFull = KINSetOptions('Verbose',true, ...
                        'MaxNewtonStep',mxnewt, ...
                        'LinearSolver', 'Dense');
        
        
    end

    % Run KINSol
    KINInit(f, neq, optionsFull);
    [status, yDash] = KINSol(yDash, strategy, yscale, fscale);
        
    % Retrive KINSol metrics
    JFNKStats = KINGetStats();
    totalJFNKLinBandTime = toc(JFNKLinBandTime);
    KINFree;

    clear optionsJNFKBandLinAlg;
    
    disp([num2str(N),'x',num2str(M),' Delta x=',...
        num2str(deltaX),', y=',num2str(deltaY),', x_0=',num2str(x0),...
        ', epsilon=',num2str(eps(1)),', F=',num2str(Fr),', PreType=',preconType,' imxy=',intMethod.x(1),'/',intMethod.y(1),' n=',num2str(n)]);
    nz = (M*(2*band+1)-band*(band+1))*(N+1)^2;
    percent = nz/(M*(N+1))^2*100;
    totalTime = totalJFNKLinBandTime + LinPreRunTime;
    nonLinStep = JFNKStats.nni;
    
    % Display the runtime stats of the JFNK method, in the form (&
    % spereator)
    % Bandwidth & Percentage of matrix full & # function evaluations &
    % # nonlinear steps & average Krylov subspace size & Preconditioner setup
    % time & Preconditioner solve time & Nonlinear function runtime & Total
    % time
    try
    avgSspace = JFNKStats.LSInfo.nli/nonLinStep;
    disp([num2str(band),' & ',num2str(percent),' & ',num2str(NonLinFuncRunCount),...
        ' & ',num2str(nonLinStep),' & ',num2str(avgSspace),...
        ' & ',num2str(LinPreRunTime),' & ',num2str(PreSolveTime),...
        ' & ',num2str(NonLinFuncRunTime),' & ',num2str(totalTime),'\\']);
    catch
        disp([num2str(band),' & ',num2str(percent),' & ',num2str(NonLinFuncRunCount),...
        ' & ',num2str(nonLinStep),' & ','N/A',...
        ' & ',num2str(LinPreRunTime),' & ',num2str(PreSolveTime),...
        ' & ',num2str(NonLinFuncRunTime),' & ',num2str(totalTime),'\\']);
    end
    
    % Save the results
    timestamp = datestr(now,'yyyy-mm-dd (HH.MM.SS)');
    
    if singType
        singularity = 'Source';
        singPara = 'eps';
    else
        singularity = 'Dipole';
        singPara = 'mu';
    end
    
    filepath = [folder,singularity,' ',num2str(N),'x',num2str(M),...
        ' ',singPara,' ',num2str(eps(1)),' F ',num2str(Fr),...
        ' x0 ',num2str(x0),' dx ',num2str(deltaX),' dy ',num2str(deltaY),...
        ' band ',num2str(band),' im-',intMethod.x(1),intMethod.y(1),' n ',num2str(n),' - ',timestamp,'.mat'];
    
    try
        save(filepath,'M','N','deltaX','deltaY','x0','source','eps','Fr','intMethod','singType','n','band','preconType','yDash0',...
        'JFNKStats','yDash','LinPreRunTime','NonLinFuncRunCount','NonLinFuncRunTime','PreSolveTime',...
        'avgSspace','fscale','maxl','mxnewt','neq','nonLinStep','nz','percent','status','strategy',...
        'timestamp','totalJFNKLinBandTime','totalTime','yscale');
    catch
        save(filepath,'M','N','deltaX','deltaY','x0','source','eps','Fr','intMethod','singType','n','band','preconType','yDash0',...
        'JFNKStats','yDash','LinPreRunTime','NonLinFuncRunCount','NonLinFuncRunTime','PreSolveTime',...
        'fscale','maxl','mxnewt','neq','nonLinStep','nz','percent','status','strategy',...
        'timestamp','totalJFNKLinBandTime','totalTime','yscale');
    end
end
end