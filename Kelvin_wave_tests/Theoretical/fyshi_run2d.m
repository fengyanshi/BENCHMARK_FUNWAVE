clear all
% Parameters------------------------------
Fr = 0.6;               % The Froude number
epsilon = 1;            % The pressure strength
freq = 10;              % The sampling frequency
deltaX = 1/freq;        % The mesh spacing
x0 = 0;                 
xEnd = 150;
deltaY = 0.1;

x = (x0:deltaX:xEnd);

% The minimum distence to the sensor (must be relatively large to ensure 
% useage of only the single integral term of equation (6.6) is valid,
% y>2/Fr^2 should suffice).
y = (0:deltaY:50);

% Calculate solutions
[zeta]=linearZetaWaveTrain(x,y,epsilon,Fr);

zeta_rotate=zeta';
[n m]=size(zeta_rotate);
n1=2*(n-1)+1;
m1=m;
zeta1=zeros(n1,m1);
zeta1(1:n,1:m)=flipud(zeta_rotate);
zeta1(n+1:n1,1:m)=zeta_rotate(2:n,1:m);

x1=[0:m1-1]*deltaX;
y1=[0:n1-1]*deltaY-(n-1)*deltaY;


figure(1)
clf
colormap jet
pcolor(x1,y1,zeta1),shading flat
colorbar
mz=max(max(zeta));
caxis([-mz mz])
xlabel('x')
ylabel('y')
name1=['Fr = ' num2str(Fr)];
eval('mkdir plots');
fname1=['plots/' 'Fr_' num2str(Fr) '.jpg'];
title(name1)
print('-djpeg',fname1)


