clear all
cases='F06_3layers';
fdir = ['/Volumes/Solid/Kelvin_wave/F2WAVE/' cases '/'];

files=[1:];

dx=0.5;
dy=0.5;

wid=8;
len=3;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
colormap jet
clf

for k=1:length(files)
fnum=sprintf('%.4d',files(k));
eta=load([fdir 'eta_' fnum],'-ASCII');

[n m]=size(eta);
x=[0:m-1]*dx;
y=[0:n-1]*dy;


pcolor(x,y,eta),shading interp
caxis([-0.5 0.5])
xlabel('x (m)')
ylabel('y (m)')

eval(['mkdir ' 'plots/' cases]);

fname=['plots/' cases '/' 'snap_' cases '_' num2str(files(k)) '.jpg'];
print('-djpeg100', fname);

end


