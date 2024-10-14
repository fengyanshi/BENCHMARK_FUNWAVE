clear all
cases='F06_L20';
fdir = ['/Volumes/Solid/Kelvin_wave/FUNWAVE/' cases '/'];

files=[1:74];

m=3072;
n=512;
dx=0.5;
dy=0.5;
x=[0:m-1]*dx;
y=[0:n-1]*dy;
[X Y]=meshgrid(x,y);

sta=load('station.txt');
xs=sta(:,1)*dx;
ys=sta(:,2)*dy;

fname=['../plots/' 'eta_' cases '.mp4' ];

% define movie file and parameters
myVideo = VideoWriter(fname,'MPEG-4');
myVideo.FrameRate = 10;  
myVideo.Quality = 100;
%vidHeight = 576; %this is the value in which it should reproduce
%vidWidth = 1024; %this is the value in which it should reproduce
open(myVideo);

fig=figure(1);
wid=10;
len=5;
set(fig,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
colormap jet


for k=1:length(files)

fnum=sprintf('%.5d',files(k));
eta=load([fdir 'eta_' fnum],'-ASCII');


clf
pcolor(X,Y,eta),shading flat
colorbar
caxis([-0.5 0.5])
xlabel('x(m)')
ylabel('y(m)')
title(['Case: ' cases '  Time = '  num2str(files(k)*2.0) ' sec'])

hold on
plot(xs,ys,'ro')
for kk=1:length(xs)
text(xs(kk)+10,ys(kk),num2str(kk))
end

pause(0.1)


% save image
F = print('-RGBImage','-r300');
%J = imresize(F,[vidHeight vidWidth]);
mov(k).cdata = F;

writeVideo(myVideo,mov(k).cdata);


end
close(myVideo)



%eval(['mkdir ../plots/' cases])
%fname=['../plots/' cases '/eta_' cases '_num_' num2str(files(k)) '.jpg' ];
%print('-djpeg100', fname);

