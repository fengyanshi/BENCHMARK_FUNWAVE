clear all

%folder_result='/Users/fengyanshi/tmp3/';
folder_result='../results_pl_shock/';
data=load('ht_setup_031041.txt');
x_data=data(:,1);
ht_data=data(:,2);
set_data=data(:,3);

nsta=119;
dx=0.025;
noise_width=5;

icount=0;
for num=1:nsta
fnum=sprintf('%.4d',num);
time_eta=load([folder_result 'sta_' fnum]);
station(:,num)=time_eta(300:end,2);
end

[ntime,nsta]=size(station);

for i=1:nsta
sta=station(:,i);
indtmp=crossing(sta);

% remove noise
ic=0;
clear ind1
for jj=2:length(indtmp)
if indtmp(jj)-indtmp(jj-1) > noise_width
ic=ic+1;
ind1(ic)=indtmp(jj);
end
end

if length(ind1)>3
jcount=0;
jmax=0;
jmin=0;
javg=0;
for j=1:2:length(ind1)-2
vcount=0.0;
vmax=0.0;
vmin=0.0;
ncount=0.0;
for jj=ind1(j):ind1(j+2)-1
ncount=ncount+1;
vcount=vcount+sta(jj);
if sta(jj)>vmax
vmax=sta(jj);
end
if sta(jj)<vmin
vmin=sta(jj);
end
end %  end jj
jcount=jcount+1;
jmax=jmax+vmax;
jmin=jmin+vmin;
vavg=vcount/ncount;
javg=javg+vavg;
end % end j
height(i)=jmax/jcount-jmin/jcount;
zavg(i)=javg/jcount;
end % endif

x_num(i)=(i-1)*dx*4.0;

end % end i


figure
subplot(211)
plot(x_data,ht_data,'ro')
hold on
plot(x_num,height,'-x')
grid
ylabel('wave height(m)')
xlabel('x (m)')
legend('data','model','Location','NorthWest')

subplot(212)
plot(x_data,set_data,'rx')
hold on
plot(x_num,zavg,'-x')
grid
ylabel('setup(m)')
xlabel('x (m)')

output(:,1)=x_num;
output(:,2)=height;
output(:,3)=zavg;

save -ASCII shock_cap_case.txt output
 
