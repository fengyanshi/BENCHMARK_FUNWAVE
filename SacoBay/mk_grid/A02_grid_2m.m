clear all

dep=load('depth_17m_dx10.txt');
dx=1.0;
dy=1.0;

[n m]=size(dep);
x=[0:m-1]*dx;
y=[0:n-1]*dy;
[X,Y]=meshgrid(x,y);

figure
clf
pcolor(X,Y,-dep),shading flat
colormap jet
colorbar
caxis([-20 20])

% coarse grid 2.0m

dx1=2.0;
dy1=2.0;

x1=[1:dx1:3500-1];
y1=[150:dy1:4700];
[X1,Y1]=meshgrid(x1,y1);
dep1=griddata(X,Y,dep,X1,Y1);

save -ASCII depth_17m_dx20.txt

% periodic ------
[n1 m1]=size(dep1);
dep4=dep1;
ntr=300;
for j=1:ntr
    for i=1:m1
        dep4(j,i)=dep1(end,i)+(dep1(ntr,i)-dep1(end,i))*(j-1)/(ntr-1.0);
    end
end

figure(3)
clf

pcolor(X1,Y1,-dep4),shading flat
colormap jet
%caxis([-15 15]);

print -djpeg100 grid_2m_1750x2276.jpg

save -ASCII depth_2m_1750x2276.txt dep4



