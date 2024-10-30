fdir = '/Users/fengyanshi/TMP/tmp2/';

m=256;
n=64;
l=3;

dx=2.0;
dy=2.0;
x=[0:m-1]*dx;
y=[0:n-1]*dy;

ns=input('ns=');
ne=input('ne=');

wid=8;
len=3;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);

colormap jet
icount=0;
for num=ns:1:ne

icount=icount+1;

fnum=sprintf('%.4d',num);
eta=load([fdir 'eta_' fnum]);
clf
pcolor(x,y,eta),shading interp
caxis([-0.5 0.5])
xlabel('x (m)')
ylabel('y (m)')
pause(1)
end
print -djpeg100 ship_wake.jpg


