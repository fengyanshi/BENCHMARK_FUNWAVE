clear all

% theory
g=9.81;
h=[1:0.5:18];
lambda=20.0;
k=2.0*pi/lambda;
kh=k.*h;
sigma=sqrt(g*k*tanh(kh));
T=2.*pi./sigma;

% numerical
data=load('../results/periods.txt');
h_num=data(:,1);
k_num=2.0*pi/lambda;
kh_num=k_num.*h_num;
T_num=data(:,4);

fig=figure(1);
clf
wid=6;
len=4;
set(fig,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
plot(kh,T,'k-',kh_num,T_num,'r--','LineWidth',2)
grid
legend('theory','FUNWAVE')
xlabel('kh')
ylabel('T(s)')

eval(['mkdir ' 'plots'])

fname=['plots/T_kh.jpg'];
print('-djpeg',fname)






