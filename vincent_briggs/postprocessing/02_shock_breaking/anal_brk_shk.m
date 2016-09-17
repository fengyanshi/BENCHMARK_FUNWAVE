clear all; clc;

fdir_res='../result_02/';
fdir_work='../02_Shock_breaking/';
%%%%%%%%%%%%
outdt=0.05;
mx=760;    ny=501;
dx=0.05;    dy=0.05;

tstart=24.9;
%tend=tstart+36.4;
tend=tstart+28.6;
H0=0.135;
xpos=27.2;
%%%%%%%%%%%%

fn=load('fname.in');
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
     if mod(kk,10) == 0
     txt=['read completed ...' num2str(floor((kk-ind1)/(ind2-ind1)*100)) '%'];
     disp(txt);
    end
    fname = ['eta_',num2str(fn(kk),'%05d'),'.mat'];
    load( [fdir_res fname ] );
    tote(:,nn)=eta(:);
    clear eta
    nn=nn+1;    
    
end
clear kk

totgrid=mx*ny;
hrms=zeros( totgrid, 1 );
for kk=1:totgrid;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pname=['result_shk',num2str(tstart),'.mat'];
save(pname);  clear pname
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hrms=reshape( hrms,[ny,mx] );
hrms=hrms/H0;

x = 0 : dx : (mx-1)*dx;
y = 0 : dy : (ny-1)*dy;
x=x';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indx=find( min(abs(x-xpos)) == abs(x-xpos) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hrmsx=hrms( :, indx );
set(figure(1),'position',[50 50 580 250])
plot( y,hrmsx,'b','LineWidth',1.4 )
xlim([9.2 15.8])
ylim([0 2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obs=load('Br_Obs.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
plot( obs(:,1), obs(:,2),'ko','Markersize',8,'Linewidth',2.0 )

xlabel('Y (m)','fontsize',12);
ylabel('H_{rms}/H_{wavemaker}','fontsize',12);

fname=['result_shk',num2str(tstart),'_',num2str(tend),'.eps'];
set(gcf,'PaperPositionMode', 'auto');
print('-depsc2',fname);

pname = ['result_shk',num2str(tstart),'_',num2str(tend),'.png'];
set(gcf,'PaperPositionMode', 'auto');
print ('-dpng', pname);
