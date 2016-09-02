clear all
clf
fdir1='/Users/fengyanshi/tmp3/';
fdir2='/Users/fengyanshi/tmp4/';

dep=load([fdir1 'dep.out']);
m=length(dep);
dx=0.025;
x=(0:m-1)*dx;

x1=30;
x2=45;

ibeg=input('ibeg');
iend=input('iend');

for num=ibeg:iend
fnum=sprintf('%.4d',num);
eta1=load([fdir1 'eta_' fnum]);
eta2=load([fdir2 'eta_' fnum]);
%nu=load([fdir1 'nubrk_' fnum]);
nu3=load([fdir2 'nubrk_' fnum]);
clf
subplot(211)
plot(x,eta1(2,:),'b')
hold on
plot(x,eta2(2,:),'r--')
plot(x,-dep(2,:),'k')
grid
axis([x1 x2 -0.1 0.15])

subplot(212)
%plot(x,nu(2,:),'k--')
hold on
plot(x,nu3(2,:),'r')
grid
axis([x1 x2 0.0 0.04])
pause

end