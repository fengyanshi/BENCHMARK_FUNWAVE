function [pAng,gof]=getPeakAngle(yDash,M,N,deltaX,deltaY,x0,Fr,varargin)
% getPeakAngle 
% yDash - vector of zetaX
% M,N - the number of nodes in the y and x directions respectively
% deltaX,deltaY - the distance between nodes in the x and y directions
% x0 - the smallest x value
% Fr - The Froude number
% varargin - (optional) is any extra input is given a contour plot of the
%            surface with the highest peaks and wake angle marked

% Check if there are additional inputs, and if so confirm plot
plotCheck = false;
if ~isempty(varargin)
    plotCheck = true;
end

% Set up surface
x = linspace(x0,x0+deltaX*(N-1),N)'/Fr^2;
y = linspace(0,deltaY*(M-1),M)/Fr^2;
[X,Y] = meshgrid(x,y);
X=X';
Y=Y';
[zeta,~,~,~] = getValues(yDash,M,N,deltaX,deltaY,x0);

% Find the location of the first downstream peak on the centreline
clinePeakpos = x([false;zeta(2:end-1,1)>zeta(1:end-2,1)&zeta(2:end-1,1)>zeta(3:end,1);false]);
clinePeakpos(clinePeakpos<0) = [];
start = clinePeakpos(1);

wavelength = 2*pi;

% Starting at the first centreline peak, check every wavelength for the
% wave peak
peakPos = [];
i=1;
for peaksX = start:wavelength:max(x);
    
    % Isolate a single transverse wave using the simplifying assumption 
    % that a transverse wave crest follow the circumference of a circle
    % whose radius is the distance from the x=0 to the centreline crest and
    % is centred two radial distances from x=0
    test = (X-2*peaksX).^2+Y.^2-peaksX^2;
    test = ~(test<=(pi^2+2*pi*peaksX)&test>=(pi^2-2*pi*peaksX)&X<(2*peaksX));
    zetaTemp = zeta;
    zetaTemp(test) = 0;
    
    % Find and store the posistion of the highest peak
    [maxZeta,inds] = max(zetaTemp);
    [~,yInd] = max(maxZeta);
    xInd = inds(yInd);
    peakPos(i,:) = [x(xInd)*Fr^2,y(yInd)*Fr^2];
    i=i+1;
end

% Discard the first and last couple of peaks when determining the wake
% angle to avoid truncation error. (peak range requires tuning)
mpeakPos = peakPos(2:(end-1),:);
[fit1,gof] = fit(mpeakPos(:,1),mpeakPos(:,2),'poly1');

% Get the angle in degrees
pAng = atan(fit1.p1)*180/pi;

% If a plot was requested
if plotCheck
    
    % Generate the surface mesh
    x = linspace(x0,x0+deltaX*(N-1),N)';
    y = linspace(0,deltaY*(M-1),M);
    [X,Y] = meshgrid(x,y);
    zeta = zeta';
      
    newzeta = [zeta(M:-1:1,:);zeta(2:M,:)];
    newY = [-Y(M:-1:1,:);Y(2:M,:)];
    newX = [X(M:-1:1,:);X(2:M,:)];
    
    % Get the edge point of the line of best fit
    z1 = fzero(fit1,0);
    y1 = feval(fit1,x(end));

    % Plot positive contours of the surface
    figure;
    maxZ = max(max(newzeta));
    contour(newX,newY,newzeta,linspace(0,maxZ,22));

    % Mark centreline, line of best fir (wake angle) and highest peaks
    hold on;
    plot([x0,x(end)],[0,0],'k--');
    plot([z1,x(end)],[0,y1],'k','LineWidth',2);
    plot(mpeakPos(:,1),mpeakPos(:,2),'r*');
    
end