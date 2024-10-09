function []=STFT_part(in1,STFT,anal_win, synth_win, hop, nfft, fs,f0,t0,t, x,index)
STFT10=STFT(:);
ci1=find(in1==0);
STFT10(ci1)=-0+0*i;
STFT1=reshape(STFT10,length(f0),length(t0));
[x_istft1, t_istft1] = istft(STFT1, anal_win, synth_win, hop, nfft, fs);
% figure
if index==3
else
plot(t, x,'color',[.5 .5 .5],'LineWidth',2)
end
grid on
xlim([0 max(t)])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
% xlabel('Time, s')
% ylabel('Signal amplitude')
% title('Original and reconstructed signal')

% plot the resynthesized signal 
hold on
if index==1
plot(t_istft1, x_istft1, '-r','LineWidth',2)
set(gca,'xticklabel','');
elseif index==2
plot(t_istft1, x_istft1, '-m','LineWidth',2)%;grid on
set(gca,'xticklabel','');
elseif index==3
plot(t_istft1, x_istft1, '-g','LineWidth',2)
ylim([-0.02 0.02])
%set(gca,'xticklabel','');
elseif index==4
plot(t_istft1, x_istft1, '-k','LineWidth',2)
set(gca,'xticklabel','');
elseif index==5
plot(t_istft1, x_istft1, '-b','LineWidth',2)

end
% legend('Original signal', 'Reconstructed signal')
