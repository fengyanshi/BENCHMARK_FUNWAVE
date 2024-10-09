function [Func,Flag] = BIFunction(variables,source,epsilon,singType,M,N,deltaX,deltaY,x0,Fr,intMethod,n)
% BIFunction The function that needs to be minimised. Takes
% the following inputs
% variables - vector of zetaX
% source - a vector of source/sink position
% epsilon - strength of the source/sink
% singType - (bool) The type of singularity, true=source, false=dipole
% M,N - the number of nodes in the y and x directions respectively
% deltaX,deltaY - the distance between nodes in the x and y directions
% x0 - the smallest x value
% Fr - the Froude number
% intMethod - struct with properties x & y storing the methods of
%             integration in their respective directions
% n - parameter for the radiation condition

% Define global variables for metrics and the GPU kernel
global NonLinFuncRunCount NonLinFuncRunTime kern gpu;

NonLinFuncRunCount = NonLinFuncRunCount + 1;
NonLinFuncRunTimeTemp = tic;

% Set success flag
Flag = 0;


% Initialise variables
Func = zeros(2*M*(N+1),1);

x = linspace(x0,x0+deltaX*(N-1),N)';
y = linspace(0,deltaY*(M-1),M);

% Use auxiliary function to take variables and return the mesh values 
[zeta,zetaX,zetaY,phi,phiX,phiY] = getValues(variables,M,N,deltaX,deltaY,x0);

% Calcualte half mesh points
xHalf = (x(1:N-1)+x(2:N))./2;
yHalf = y;

zetaHalf = (zeta(1:N-1,:)+zeta(2:N,:))./2;
zetaXHalf = (zetaX(1:N-1,:)+zetaX(2:N,:))./2;
zetaYHalf = (zetaY(1:N-1,:)+zetaY(2:N,:))./2;
phiHalf = (phi(1:N-1,:)+phi(2:N,:))./2;
phiXHalf = (phiX(1:N-1,:)+phiX(2:N,:))./2;
phiYHalf = (phiY(1:N-1,:)+phiY(2:N,:))./2;

% Use auxiliary function to compute correct integration weightings
weightingX=intWeight(x,intMethod.x);
weightingY=intWeight(y,intMethod.y);

if gpu
% Compute function on the GPU
dFunc = feval(kern,N*ones(1,1,'int64'),M*ones(1,1,'int64'),...
    zeta,zetaX,zetaY,zetaHalf,zetaXHalf,zetaYHalf,...
    phi,phiX,phiHalf,phiXHalf,phiYHalf,...
    x,y',xHalf,yHalf',weightingX,weightingY',...
    epsilon,source',singType,Fr,...
    n,Func);
% Copy results back to host memory
Func = gather(dFunc);
else
    weiX = repmat(weightingX,1,M);
    % For every half point in the mesh
    % Free Surface condition
    Func1=1/2*((1+zetaXHalf.^2).*phiYHalf.^2+(1+zetaYHalf.^2).*phiXHalf.^2-2*zetaXHalf.*zetaYHalf.*phiXHalf.*phiYHalf)./(1+zetaXHalf.^2+zetaYHalf.^2)+zetaHalf/Fr^2-1/2;
    Func1 = Func1(:);
    Func2 = zeros(M*(N-1),1);
    for i=1:(M*(N-1))
            l = floor((i-1)/(N-1))+1;
            k = mod(i-1,N-1)+1;
            % Initialise sums
            A = 1 + zetaXHalf(k,l)^2;
            B = 2*zetaXHalf(k,l)*zetaYHalf(k,l);
            C = 1 + zetaYHalf(k,l)^2;

                % Calculate often used values
                xDiff = (x-xHalf(k)*ones(N,1));
                xDiffSq = xDiff.^2;
                xDiff = repmat(xDiff,1,M);
                xDiffSq = repmat(xDiffSq,1,M);

                yNegDiff = (y-yHalf(l)*ones(1,M));
                yNegDiffSq = yNegDiff.^2;
                yNegDiff = repmat(yNegDiff,N,1);
                yNegDiffSq = repmat(yNegDiffSq,N,1);

                yPosDiff = (y+yHalf(l)*ones(1,M));
                yPosDiffSq = yPosDiff.^2;
                yPosDiff = repmat(yPosDiff,N,1);
                yPosDiffSq = repmat(yPosDiffSq,N,1);

                zetaDiff = (zeta-zetaHalf(k,l)*ones(N,M));
                radiusSquaredYNeg = sqrt(xDiffSq+yNegDiffSq+zetaDiff.^2);
                radiusSquaredYPos = sqrt(xDiffSq+yPosDiffSq+zetaDiff.^2);

                % Calculate complicated values for the integral
                K1 = (zetaDiff-xDiff.*zetaX-yNegDiff.*zetaY)./...
                        (radiusSquaredYNeg).^3+...
                        (zetaDiff-xDiff.*zetaX-yPosDiff.*zetaY)./...
                        (radiusSquaredYPos).^3;
                K2 = 1./radiusSquaredYNeg+1./radiusSquaredYPos;
                S2 = 1./sqrt(A*xDiffSq+B*xDiff.*yNegDiff+C*yNegDiffSq)+...
                        1./sqrt(A*xDiffSq-B*xDiff.*yPosDiff+C*yPosDiffSq);

                xs = repmat(x,1,M);
                I1all = (phi-xs-(phiHalf(k,l)-xHalf(k))*ones(N,M)).*K1;
                I2all = (zetaX.*K2-zetaXHalf(k,l)*S2);

                % Sum over all x
                I1Inner = sum(weiX.*I1all);
                I2InnerDash = sum(weiX.*I2all);

                % Sum over all y
                I1 = sum(weightingY.*I1Inner);
                I2Dash = sum(weightingY.*I2InnerDash);

            % Calculate often used values
            sUpper = x(N)-xHalf(k);
            sLower = x(1)-xHalf(k);
            tNegUpper = y(M)-yHalf(l);
            tNegLower = y(1)-yHalf(l);
            tPosUpper = -y(M)-yHalf(l);
            tPosLower = -y(1)-yHalf(l);

            % Calculate parts of the closed form integral
            singInt1 = definitInt(sLower,sUpper,tNegLower,tNegUpper,A,B,C);
            singInt2 = definitInt(sLower,sUpper,tPosLower,tPosUpper,A,B,C);

            % Join parts
            singInt = singInt1 - singInt2;

            I2DoubleDash = zetaXHalf(k,l)*singInt;

            % Get the contribution from the singularity
            if singType
                singCont = epsilon./sqrt((xHalf(k)-source(1)).^2+...
                    (yHalf(l)-source(2)).^2+...
                    (zetaHalf(k,l)-source(3)).^2);
            else
                singCont = -epsilon*(xHalf(k)-source(1))./sqrt((xHalf(k)-source(1)).^2+...
                    (yHalf(l)-source(2)).^2+...
                    (zetaHalf(k,l)-source(3)).^2).^3;
            end

            % Calculate the value of the function
            Func2(i)= 2*pi*(phiHalf(k,l)-xHalf(k)) +...
                singCont - I1 - I2Dash - I2DoubleDash;

    end
    % Apply radiation conditions
    Func3 = x0*(phiX(1,:)-1)+n*(phi(1,:)-x0);
    Func4 = x0/deltaX*(phiX(2,:)-phiX(1,:))+n*(phiX(1,:)-1);
    Func5 = x0*zetaX(1,:)+n*zeta(1,:);
    Func6 = x0/deltaX*(zetaX(2,:)-zetaX(1,:))+n*zetaX(1,:);
    
    % Arrange equations in the output
    for i=1:M
        pos = (i-1)*(N-1)+1;
        pos1 = (i-1)*(N+1)+1;
        pos2 = (i-1)*(N+1)+1+M*(N+1);
        Func(pos1) = Func3(i);
        Func(pos1+1) = Func4(i);
        Func((pos1+2):(pos1+2+N-2))=Func1(pos:(pos+N-2));
        Func(pos2) = Func5(i);
        Func(pos2+1) = Func6(i);
        Func((pos2+2):(pos2+2+N-2))=Func2(pos:(pos+N-2));
    end

end

NonLinFuncRunTime = NonLinFuncRunTime + toc(NonLinFuncRunTimeTemp);
end

% Determines one half of the closed form integral
function val=closedIntt(s,t,A,B,C)

if abs(t)<10^-15
    val = 0;
else
    val= t/sqrt(A)*log(2*A*s+B*t+2*sqrt(A*(A*s^2+B*s*t+C*t^2)));
end

end

% Determines another half of the closed form integral
function val=closedInts(s,t,A,B,C)

if abs(s)<10^-15
    val = 0;
else
    val= s/sqrt(C)*log(2*C*t+B*s+2*sqrt(C*(A*s^2+B*s*t+C*t^2)));
end

end

% Determines the closed form integral
function singInt=definitInt(sLower,sUpper,tLower,tUpper,A,B,C)

singInt =   (closedIntt(sUpper,tUpper,A,B,C)+closedInts(sUpper,tUpper,A,B,C))-...
            (closedIntt(sLower,tUpper,A,B,C)+closedInts(sLower,tUpper,A,B,C))-...
            (closedIntt(sUpper,tLower,A,B,C)+closedInts(sUpper,tLower,A,B,C))+...
            (closedIntt(sLower,tLower,A,B,C)+closedInts(sLower,tLower,A,B,C));

end

