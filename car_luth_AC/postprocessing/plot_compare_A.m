clear all; clc; clf;

tskp=45.15;
%time=load('..\time_dt.out');
%time=time(:,1);

set(gcf,'units','inches','paperunits','inches','papersize', [6 8],'position',[1 1 6 8],'paperposition',[0 0 6 8]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lab={'22' '24' '30.5' '32.5' '33.5' '34.5' '35.7' '37.3' '39.0' '41'};

for num_fig=1:10
    subplot(5,2,num_fig)
    file_number = sprintf('%.4d',num_fig);
    eta=load(['../01_CaseA/sta_' file_number]);
    data=load(['../lab_data/CaseA/Exp_' lab{num_fig} 'm.txt']);
    
    
    plot( eta(:,1),eta(:,2), 'k-', 'Linewidth', 1.5 );  hold on;
plot( data(:,1)+tskp, data(:,2), 'ko', 'MarkerSize', 4., 'LineWidth',1.2 )

ytickname=[-0.04 -0.02 0.0 0.02 0.04];
set(gca,'ytick',ytickname,'fontsize',11);
grid

ylim([ -0.03 0.035 ])
xlim([ 45 49.8])

name=['x= ' lab{num_fig} '(m)'];
    title(name)

ylabel(['\eta (m)'],'fontsize',10)
set( gca, 'fontsize', 10 )
if num_fig >=9    
xlabel('time (s)')
end
    clear eta data;
    
    
end

print -djpeg model_data_A.jpg
break

