clear all
cases='F12_3layers';
fdir = ['/Volumes/Solid/Kelvin_wave/F2WAVE/' cases '/'];


site=[1:16];

df=0.0005;

for k=1:length(site)
xblock1{site(k)}=[20,60,60,20,20];yblock1{site(k)}=[0.12,0.22,0.34-df,0.22-df,0.14];
xblock2{site(k)}=[35,53,53,35,35];yblock2{site(k)}=[0.27,0.325,0.40,0.36,0.27];
xblock3{site(k)}=[550,980,980,550,550];yblock3{site(k)}=[0.0025,0.0025,0.0075,0.0075,0.0025];
xblock4{site(k)}=[200,980,980,200,200];yblock4{site(k)}=[0.0075+df,0.0075+df,0.0125-df,0.0125-df,0.0075+df];
xblock5{site(k)}=[200,980,980,200,200];yblock5{site(k)}=[0.0125,0.0125,0.02,0.02,0.0125];
end

for k=1:length(site)
fnum=sprintf('%.4d',site(k));
data=load([fdir 'sta_' fnum],'-ASCII');

time=data(1:1:600,1);  % use 625 for f12. 490 for f14, 395 for f16
eta=data(1:1:600,2);

dt=0.1;
fs=1/dt;

end_time=time(end);
new_time=[10*dt:dt:end_time];

x=interp1(time,eta,new_time);

xlen = length(x);                   
t = (0:xlen-1)/fs;                  
% Input:
% x - signal in the time domain
% win - analysis window function
% hop - hop size
% nfft - number of FFT points
% fs - sampling frequency, Hz

% define the analysis and synthesis parameters
wlen = 512;
hop = wlen/32;
nfft = 4*wlen;

% generate analysis and synthesis windows
anal_win = blackmanharris(wlen, 'periodic');
synth_win = hamming(wlen, 'periodic');

% perform time-frequency analysis and resynthesis of the signal
[STFT,f0,t0] = stft(x, anal_win, hop, nfft, fs);
[x_istft, t_istft] = istft(STFT, anal_win, synth_win, hop, nfft, fs);
psd = abs(STFT).^2/xlen;
psd = psd./fs;


len_win=400;
[s9,f9,t9,ps9,fc9,tc9]=spectrogram(x,len_win,len_win-10,len_win,fs);%,length(xs(:,1)),'power','yaxis');
figure(2);
clf
set(gcf,'Color','w','Position', [950 550 500 650]);
% set(gca,'LooseInset',get(gca,'TightInset'));
pos1 = [0.18 0.68 0.675 0.30];
subplot('Position',pos1)
hold on
% [c,h]=contourf(t0,f0,log10(abs((STFT))));
[c,h]=contourf(t9,f9,log10(ps9),[-10:0.1:3]);
%set(h,'LineColor','none')
ch=colorbar('position',[0.87 0.68 0.02 0.3]);
set(get(ch,'title'),'string','$$log10(PSD)$$','interpreter','latex','Rotation',90.0,'FontSize',12);%,'position',[590 -15]);
caxis([-10,3]);
lbpos = get(ch,'title');
set(lbpos, 'position', [50,80,0]) ;

%set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
%set(gca,'xticklabel','');
ylabel('$$f(Hz)$$','interpreter','latex','FontSize',12)
%%
% ---------
xv1=xblock1{site(k)}; yv1=yblock1{site(k)};
xv2=xblock2{site(k)}; yv2=yblock2{site(k)};
xv3=xblock3{site(k)}; yv3=yblock3{site(k)};
xv4=xblock4{site(k)}; yv4=yblock4{site(k)};
xv5=xblock5{site(k)}; yv5=yblock5{site(k)};

hold on
plot(xv1,yv1,'r-','linewidth',2);
plot(xv2,yv2,'g-','linewidth',2);
%plot(xv3,yv3,'m-','linewidth',2);
%plot(xv4,yv4,'k-','linewidth',2);
%plot(xv5,yv5,'b-','linewidth',2);
[X,Y] = meshgrid(t0,f0);
xq=X(:);yq=Y(:);
[in1,on1] = inpolygon(xq,yq,xv1,yv1);
[in2,on2] = inpolygon(xq,yq,xv2,yv2);
[in3,on3] = inpolygon(xq,yq,xv3,yv3);
[in4,on4] = inpolygon(xq,yq,xv4,yv4);
[in5,on5] = inpolygon(xq,yq,xv5,yv5);
% % plot(xq(in1),yq(in1),'ro') % points inside
% % plot(xq(in2),yq(in2),'g+') % points inside
% % plot(xq(in3),yq(in3),'y+') % points inside
%%
%axis([0 max(t) 0 1.0]);box on
plot([0 100],[0.3 0.3],'w--','LineWidth',2)
text(22,0.35,'kh = \pi','Color','w','FontSize',13)

ax1=[20 180 0 0.6];
ax2=[20 180 -0.5 0.5];

axis(ax1);box on

% part ----------
pos2 = [0.18 0.40 0.675 0.2];
pos4 = [0.18 0.15 0.675 0.2];
subplot('Position',pos2)
STFT_part(in1,STFT,anal_win, synth_win, hop, nfft, fs,f0,t0,t, x,1);
 ylabel('$$\eta(m)$$','interpreter','latex','FontSize',13)
legend('time series','divergence chirp')
axis(ax2);box on
%
subplot('Position',pos4)
STFT_part(in2,STFT,anal_win, synth_win, hop, nfft, fs,f0,t0,t, x,3);
ylabel('$$\eta(m)$$','interpreter','latex','FontSize',13)
legend('fake high frequencies')
axis(ax2);box on
%
%set(gca,'xtick',[])
% ylabel('$$\eta(m)$$','interpreter','latex','FontSize',13)

xlabel('$$Time(s)$$','interpreter','latex','FontSize',13);



% ----------
%title(['Case: ' cases(1:3) ' ' cases(5:end) ', gauge: ' num2str(site(k))])

fname=['../plots/' 'spc_' cases '_site_' num2str(site(k)) '.jpg' ];
print('-djpeg100', fname);

end


