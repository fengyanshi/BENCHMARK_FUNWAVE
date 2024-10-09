clear all
cases='Hres_F09_3layers';
fdir = ['/Users/fyshi/TMP/tmp2/' cases '/'];

files=[16];


for k=1:length(files)
fnum=sprintf('%.4d',files(k));
data=load([fdir 'probe_' fnum],'-ASCII');

num=1012;

time=data(1:4:num,1);
eta=data(1:4:num,2);

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

% plot the original signal
figure(1)
clf
set(gcf,'Color','w','Position', [950 550 700 750]);
set(gca,'LooseInset',get(gca,'TightInset'));
plot(t, x, 'b')
grid on
xlim([0 max(t)])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('Time, s')
ylabel('Signal amplitude')
title('Original and reconstructed signal')

% plot the resynthesized signal 
hold on
plot(t_istft, x_istft, '-.r')
legend('Original signal', 'Reconstructed signal')
%%

len_win=400;
[s9,f9,t9,ps9,fc9,tc9]=spectrogram(x,len_win,len_win-10,len_win,fs);%,length(xs(:,1)),'power','yaxis');
figure(2);
clf
set(gcf,'Color','w','Position', [950 550 700 750]);
% set(gca,'LooseInset',get(gca,'TightInset'));
pos1 = [0.18 0.65 0.675 0.30];
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
axis([0 max(t) 0 1.0]);box on
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
%set(gca,'xticklabel','');
ylabel('$$f(Hz)$$','interpreter','latex','FontSize',12)
%%

title(['Case: ' cases(1:3) ' ' cases(5:end) ', gauge: ' num2str(files(k))])

eval(['mkdir ' 'plots/' cases]);

fname=['plots/' cases '/' 'sp_' cases '_' num2str(files(k)) '.jpg'];
print('-djpeg100', fname);

end


