clear all; clc;

fdir_res='../result_03/';
fdir_work='../03_nobreaking/';
%%%%%%%%%%%%
mx=660;    ny=500;
dx=0.05;    dy=0.05;

tstart=16.9;
tend=tstart+36.4;
H0=0.0254;
xpos=22.2;
%%%%%%%%%%%%

time=load([fdir_work 'time_dt.out']);
time=time(:,1);

tmp=abs( time-tstart );
ind1=find( tmp == min(tmp) );
clear tmp

tmp=abs( time-tend );
ind2=find( tmp == min(tmp) );
clear tmp

t_intv=time(ind1:ind2);

tote=zeros( mx*ny,ind2-ind1+1 );
nn=1;
for kk=ind1:ind2
%     kk
    if mod(kk,10) == 0
     txt=['read completed ...' num2str(floor((kk-ind1)/(ind2-ind1)*100)) '%'];
     disp(txt);
    end

    fname = [ 'eta_',num2str(kk,'%05d') ];
    eta=load( [ fdir_res, fname ] );
    tote(:,nn)=eta(:);
    clear eta
    nn=nn+1;    
    
end
clear kk

totgrid=mx*ny;
hrms=zeros( totgrid, 1 );
for kk=1:totgrid;
%     kk
    if mod(kk,100) == 0
     txt=['completed ...' num2str(floor(kk/totgrid*100)) '%'];     
     disp(txt);
    end

    tmp=tote(kk,:);
    [ wh, ampc, ampt, per, nw ] = zeroup( tmp, t_intv );
    clear tmp amp* per nw
    
    h2=wh.*wh;
    sumh2=sum(h2);
    nn=length(h2);
    hrmsloc=sqrt(sumh2/nn);
    
    hrms(kk,1)=hrmsloc;    
end

clear tote
pname=['result_nobrk',num2str(tstart),'.mat'];
save(pname);  clear pname

hrms=reshape( hrms,[ny,mx] );
hrms=hrms/H0;

x = 0 : dx : (mx-1)*dx;
y = 0 : dy : (ny-1)*dy;
x=x';  %y=fliplr(y);

indx=find( min(abs(x-xpos)) == abs(x-xpos) );
hrmsx=hrms( :, indx );

set(figure(1),'position',[50 50 580 350])
plot( y,hrmsx,'k','LineWidth',1.4 )
xlim([9 16]);
ylim([0 2.5])
grid on

obs=load('Vin_Nonbr.txt');
hold on
plot( obs(:,1), obs(:,2),'ko','Markersize',8,'Linewidth',2.5 )

xlabel('Y (m)','fontsize',12);
ylabel('H_{rms}/H_{wavemaker}','fontsize',12);

% lgd=legend('FUNWAVE', ...
%     'Vincent & Briggs (1989) Observation', ...
%     'Location','NorthWest','Orientation','vertical');
% legend boxoff
% set(lgd,'fontsize',12)

fname=['result_nobrk',num2str(tstart),'_',num2str(tend),'.eps'];
set(gcf,'PaperPositionMode', 'auto');
print('-depsc2',fname);

pname = ['result_nobrk',num2str(tstart),'_',num2str(tend),'.png'];
set(gcf,'PaperPositionMode', 'auto');
print ('-dpng', pname);