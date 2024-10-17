clear all

period=[2.5 3 6];

fig=figure(1);
clf
wid=8;
len=6;
set(fig,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);

for kk=1:length(period)

cases=['T' num2str(period(kk)) 's'];
fdir = ['/Users/fyshi/OUTSIDE_Google/GITHUB_M3/FUNWAVE-TVD/simple_cases/vertical_velocity_structure/results/' cases '/'];

if period(kk)==2.5
m0=645;
end

if period(kk)==3
m0=650;
end
 
if period(kk)==6
m0=650;
end


files=[8];

if period(kk)==9
m0=1148;
files=[5];
end


dep=load([fdir 'dep.out']);
[n m]=size(dep);

z_num=10;

% theory
amp=0.5;
g=9.81;
h=10.0;
T=period(kk);
f=1/T;
om=f*2*pi;
K=wvnum_omvec(h,om,g);
lambda=2*pi/K;
C=lambda*f;
kh=K*h;
zt=-[0:h];
u_theory=amp*om*cosh(kh+K.*zt)/sinh(kh);



for k=1:length(files)

fnum=sprintf('%.5d',files(k));

eta=load([fdir 'eta_' fnum],'-ASCII');
mask=load([fdir 'mask_' fnum],'-ASCII');
u=load([fdir 'u_' fnum],'-ASCII');
v=load([fdir 'v_' fnum],'-ASCII');
Ax=load([fdir 'Ax_' fnum],'-ASCII');
Ay=load([fdir 'Ay_' fnum],'-ASCII');
Bx=load([fdir 'Bx_' fnum],'-ASCII');
By=load([fdir 'By_' fnum],'-ASCII');

% ---------------------------------
% u(z)=(za-z)Ax+0.5(za^2-z^2)Bx
% v(z)=(za-z)Ay+0.5(za^2-z^2)By
% ---------------------------------

clear z U V
for l=1:z_num
z(l,:,:)=-dep(:,:)*(l-1)/(z_num-1);
end

za=-0.5528*dep+0.4472*eta;
for l=1:z_num
zl=squeeze(z(l,:,:));
U(l,:,:)=u+(za-zl).*Ax+0.5*(za.^2-zl.^2).*Bx;
V(l,:,:)=v+(za-zl).*Ay+0.5*(za.^2-zl.^2).*By;
end

n0=25;

subplot(1,3,kk)
plot(U(:,n0,m0),z(:,n0,m0),'b-','LineWidth',2)
xlabel('u (m/s)')



grid

hold on
plot(u_theory,zt,'r--','LineWidth',2)

axis([-0.3 1.3 -10 0])

if kk==1
ylabel('z (m)')
end

if kk==1
legend('FUNWAVE','Theory','LOCATION','SouthEast')
end

title(['kh = ' num2str(kh)])

end
end

eval(['mkdir ' 'plots'])

fname=['plots/compare.jpg'];
print('-djpeg',fname)






