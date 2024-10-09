function Flag = LinearBlockJacPre(variables,Fx,M,N,deltaX,deltaY,x0,Fr,intMethod,n,yscale,fscale,band)

global LinPreRunTime data mkl;


LinAlgPreRunTimeTemp = tic;

Flag = 0;


%% Initialise variables/ generate the mesh
ps = M*(N+1);
x = linspace(x0,x0+deltaX*(N-1),N)';
y = linspace(0,deltaY*(M-1),M);


xHalf = (x(1:N-1)+x(2:N))./2;
yHalf = y;

moreX = repmat(x(2:end)',1,M);
moreY = repmat(y,N-1,1);
moreY = reshape(moreY,1,(N-1)*M);

moreXHalf = repmat(xHalf,M,1);
moreYHalf = repmat(yHalf,N-1,1);
moreYHalf = moreYHalf(:);


weightingX=intWeight(x,intMethod.x);
weiX = repmat(weightingX,1,M);
weiXless = weiX(2:end,:);
weiXless = weiXless(:);

weightingY=intWeight(y,intMethod.y);
weiY = repmat(weightingY,N,1);
weiYless = weiY(2:end,:);
weiYless = weiYless(:);

weiXless = weiXless';
weiYless = weiYless';


%% Sub matrix A
temp1 = ones(2,N-1)/2;
temp = [0,n,0;...
        1/2,n-x0/deltaX,x0;
        temp1',[x0/deltaX;zeros(N-2,1)]];
temp(N+1,1) = 0;
Asmall = spdiags(temp,[-1,0,1],N+1,N+1);
if mkl
    temp = [zeros(1,N+1);[0,0;n,n-x0/deltaX],temp1;x0,x0/deltaX,zeros(1,N-1)];
    A = repmat(temp,1,M);
else
    A = spdiags(repmat(temp,M,1),[-1,0,1],M*(N+1),M*(N+1));
end

%% Sub matrix B/C base
base = deltaX*tril(ones(N-1)-3/4*eye(N-1)-1/4*diag(ones(N-2,1),-1));
base = [zeros(2,N+1);
    ones(N-1,1),[deltaX/4;deltaX/2*ones(N-2,1)],base];

Bcoef = 1/Fr^2;
Ccoef = 2*pi;
%% Sub matrix D

D = zeros(ps,ps);
for l=1:M
    bandStart = max((l-band-1)*(N+1)+1,1);
    bandEnd = min((l+band)*(N+1),M*(N+1));
    for k = 0:(N)
        pos = (l-1)*(N+1)+k+1;
        if k==0
            D(pos,pos) = n;
        else
            KtempCol = K3(x(k)*ones(M*(N-1),1),y(l)*ones(M*(N-1),1),moreXHalf,moreYHalf);
            intCol = weightingX(k)*weightingY(l)*KtempCol;
            
            
            for i=1:M
                s1 = (i-1)*(N+1)+1+2;
                e1 = i*(N+1);
                s2 = (i-1)*(N-1)+1;
                e2 = i*(N-1);
                D(s1:e1,pos) = -intCol(s2:e2);
                if i==l
                    if k==1
                        D(s1-2,pos) = x0;
                        D(s1-1,pos) = n-x0/deltaX;
                    elseif k==2
                        D(s1-1,pos) = x0/deltaX;
                    end
                end
            end

            if k>1
                Ktemp = K3(moreX,moreY,xHalf(k-1)*ones(1,M*(N-1)),yHalf(l)*ones(1,M*(N-1)));
                Ktempextra = K3(x0*ones(1,M),y,xHalf(k-1)*ones(1,M),yHalf(l)*ones(1,M));
                int1 = weiXless.*weiYless.*Ktemp;
                int1extra = weightingX(1)*weightingY.*Ktempextra;
                S2int = sum(int1,2)+sum(int1extra,2);

                sUpper = x(N)-xHalf(k-1);
                sLower = x(1)-xHalf(k-1);
                tNegUpper = y(M)-yHalf(l);
                tNegLower = y(1)-yHalf(l);
                tPosUpper = -y(M)-yHalf(l);
                tPosLower = -y(1)-yHalf(l);

                % Calculate parts of the closed form integral
                singInt1 = definitInt(sLower,sUpper,tNegLower,tNegUpper,1,0,1);
                singInt2 = definitInt(sLower,sUpper,tPosLower,tPosUpper,1,0,1);

                S2close = singInt1 - singInt2;
                S2 = 1/2*S2int - 1/2*S2close;

                D(pos,pos) = D(pos,pos) + S2;
                D(pos,pos-1) = D(pos,pos-1) + S2;

            end
            D(1:(bandStart-1),pos)=0;
            D((bandEnd+1):end,pos)=0;
        end
    end

end

%% Form the Schur compliment
Dtemp = Ccoef*Bcoef*base*(Asmall\base);
for i=1:M
    s = (i-1)*(N+1)+1;
    e = i*(N+1);
    D(s:e,s:e) = D(s:e,s:e)-Dtemp;
end

%% Factor matrices
ipvtA = zeros(ps,1,'int64');
ipvtD = zeros(ps,1,'int64');
if mkl
factorBand(A,ipvtA,1,1);
factorDense(D,ipvtD);
end

%% Assign the proconditioner to a global variable
data.A = A;
data.ipvtA = ipvtA;
data.Bcoef = Bcoef;
data.Ccoef = Ccoef;
data.D = D;
data.ipvtD = ipvtD;

LinPreRunTime = LinPreRunTime + toc(LinAlgPreRunTimeTemp);
end

% Kernel function
function K = K3(x,y,xHalf,yHalf)

xDiff = x-xHalf;
yNegDiff = y-yHalf;
yPosDiff = y+yHalf;

K=1./sqrt(xDiff.^2+yNegDiff.^2)+1./sqrt(xDiff.^2+yPosDiff.^2);


end


% Determines one half of the closed form integral
function val=closedIntt(s,t,A,B,C)

    val= t/sqrt(A).*log(2*A*s+B*t+2*sqrt(A*(A*s.^2+B*s.*t+C*t.^2)));
    val(abs(t)<10^-15) = 0;
    
end

% Determines another half of the closed form integral
function val=closedInts(s,t,A,B,C)

    val= s/sqrt(C).*log(2*C*t+B*s+2*sqrt(C*(A*s.^2+B*s.*t+C*t.^2)));
    val(abs(s)<10^-15) = 0;
    
end

% Determines the closed form integral
function singInt=definitInt(sLower,sUpper,tLower,tUpper,A,B,C)

singInt =   (closedIntt(sUpper,tUpper,A,B,C)+closedInts(sUpper,tUpper,A,B,C))-...
            (closedIntt(sLower,tUpper,A,B,C)+closedInts(sLower,tUpper,A,B,C))-...
            (closedIntt(sUpper,tLower,A,B,C)+closedInts(sUpper,tLower,A,B,C))+...
            (closedIntt(sLower,tLower,A,B,C)+closedInts(sLower,tLower,A,B,C));

end


