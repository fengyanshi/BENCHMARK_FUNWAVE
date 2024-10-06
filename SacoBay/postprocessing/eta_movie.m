clear all
cases='result_TMA_10';
files=[359]; 
xp1=940;
xp2=930;
xp3=920;

fdir_results=['/Users/fengyanshi/OUTSIDE_Google/Users/Saco/' cases '/'];

n=2276;
m=1750;

DimsX={[m n]};

% dimensions
dx=2.0;
dy=2.0;
x=[0:m-1]*dx;
y=[0:n-1]*dy;


% define movie file and parameters
myVideo = VideoWriter('videoOut.mp4','MPEG-4');
myVideo.FrameRate = 5;  
myVideo.Quality = 100;
%vidHeight = 576; %this is the value in which it should reproduce
%vidWidth = 1024; %this is the value in which it should reproduce
open(myVideo);

fig=figure(1);
colormap jet

wid=8;
len=10;
set(fig,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);

fname=[fdir_results 'dep.out'];
fileID=fopen(fname);
dep=fread(fileID,DimsX{1},'*single');
fclose(fileID);
dep=dep';

for k=1:length(files) 

numb=files(k);

fnum=sprintf('%.5d',numb);

% read files -----------------------

fname=[fdir_results 'eta_' fnum];
fileID=fopen(fname);
eta=fread(fileID,DimsX{1},'*single');
fclose(fileID);
eta=eta';


fname=[fdir_results 'mask_' fnum];
fileID=fopen(fname);
mask=fread(fileID,DimsX{1},'*single');
fclose(fileID);
mask=mask';

eta(mask<1)=NaN;

clf


pcolor(x(1:4:end),y(1:4:end),eta(1:4:end,1:4:end)),shading flat

time1=['surface (m) at ' num2str(numb*10.0) ' s'];
title(time1)

pause(0.1)
caxis([-1.5 3.5]);

% save image
F = print('-RGBImage','-r300');
%J = imresize(F,[vidHeight vidWidth]);
mov(k).cdata = F;

writeVideo(myVideo,mov(k).cdata);

end
close(myVideo)

eval('mkdir plots');
filename=['plots/snap_' fnum];
print('-djpeg', filename)


