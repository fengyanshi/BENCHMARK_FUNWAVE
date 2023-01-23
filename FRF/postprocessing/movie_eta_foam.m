clear all

% folder of results
fdir_results='/Users/fengyanshi/TMP/tmp2/frf/';

% dimensions
m=660;
n=800;
DimsX={[m n]};

% depth
fname=[fdir_results 'dep.out'];
fileID=fopen(fname);
dep=fread(fileID,DimsX{1},'*single');
fclose(fileID);
dep=dep';
dep=flipud(dep);
dep=fliplr(dep);


% image from google earth, not the frf bathy has a rotation, cannot use google_map
RGB = imread('frf_03.jpg');

% match (x,y) and image approximately
x=[1:m]*1.+105;
y=[n:-1:1]*1.0+180;
[X Y]=meshgrid(x,y);

% define movie file and parameters
myVideo = VideoWriter('videoOut.mp4','MPEG-4');
myVideo.FrameRate = 10;  
myVideo.Quality = 100;
vidHeight = 576; %this is the value in which it should reproduce
vidWidth = 1024; %this is the value in which it should reproduce
open(myVideo);


% file range
files=[1:1:499];

fig=figure(1);
colormap jet

for k=1:length(files) 

numb=files(k);

fnum=sprintf('%.5d',numb);

% read files -----------------------
fname=[fdir_results 'eta_' fnum];
fileID=fopen(fname);
eta=fread(fileID,DimsX{1},'*single');
fclose(fileID);

eta=eta';
eta=flipud(eta);
eta=fliplr(eta);

fname=[fdir_results 'mask_' fnum];
fileID=fopen(fname);
mask=fread(fileID,DimsX{1},'*single');
fclose(fileID);
mask=mask';
mask=flipud(mask);
mask=fliplr(mask);

fname=[fdir_results 'age_' fnum];
fileID=fopen(fname);
age=fread(fileID,DimsX{1},'*single');
fclose(fileID);
age=age';
age=flipud(age);
age=fliplr(age);

% read over -------

% make up breaker and foam for visualization 

% bubbles
breaker=age;
breaker(breaker>10)=0.0;  % bubble last 10 sec

bubble=breaker;
bubble(bubble==0.0)=9000.0;
bubble=1.5*exp(-bubble/10);
bubble(bubble<0.1)=NaN;

% foam
age(age>600.0)=0.0;
age(mask<1)=NaN;
foam=age;
foam(foam==0.0)=9000.0;
foam=1.0*exp(-foam/50);

% eta without breaker
eta(mask<1)=NaN;
eta(breaker>0)=NaN;

clf

subplot(121)
B=imrotate(RGB,1);
A=imagesc(B);
hold on
pcolor(X,Y,eta),shading flat
pcolor(X,Y,bubble),shading flat

caxis([-1 3.5])
plot([110 500],[620 620],'w:','LineWidth',2)
axis([30 600 250 950])

subplot(122)
B=imrotate(RGB,1);
A=imagesc(B);
hold on
pcolor(X,Y,foam),shading flat
plot([110 500],[620 620],'w:','LineWidth',2)
axis([30 600 250 950])
caxis([0 2])

pause(0.1)


% save image
F = print('-RGBImage','-r300');
J = imresize(F,[vidHeight vidWidth]);
mov(k).cdata = J;

writeVideo(myVideo,mov(k).cdata);

end

close(myVideo)




