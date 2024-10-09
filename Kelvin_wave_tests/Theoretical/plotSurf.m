function plotSurf(yDash,M,N,deltaX,deltaY,x0)
% plotSurf Plots the surface that is the solution to the problem.
% yDash - vector of zetaX
% M,N - the number of nodes in the y and x directions respectively
% deltaX,deltaY - the distance between nodes in the x and y directions
% x0 - the smallest x value



% Initialise variables
x = linspace(x0,x0+deltaX*(N-1),N)';
y = linspace(0,deltaY*(M-1),M)';

[zeta,~,~,~,~,~] = getValues(yDash,M,N,deltaX,deltaY,x0);
% Orientate zeta matrix and create grid
zeta=zeta';
[X,Y] = meshgrid(x,y);

% Reflect values about the y axis

newzeta = [zeta(M:-1:1,:);zeta(2:M,:)];
newY = [-Y(M:-1:1,:);Y(2:M,:)];
newX = [X(M:-1:1,:);X(2:M,:)];


%% Plot Surface
figure;
surf(newX,newY,real(newzeta));
y  = [0 255];  % luminance
cb = 255*[1 1]; % chrominance-blue
cr = 0*[1 1];  % chrominance-red
map = [linspace(y(1),  y(2),  64)
       linspace(cb(1), cb(2), 64)
       linspace(cr(1), cr(2), 64)]' / 255; % in ycbcr space
rgbmap = ycbcr2rgb(map); % in rgb space

colormap(rgbmap);
caxis([-0.11 0.085]);
xlabel('x','FontName','Times New Roman','FontSize',15);
ylabel('y','FontName','Times New Roman','FontSize',15);
zlabel('\zeta','FontName','Times New Roman','FontSize',15,'Rotation',0);
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',15);
ea = min(1,70/(N-1));
c=get(gca,'Children');
set(c,'EdgeAlpha',ea);
