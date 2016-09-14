clear all

fdir1='/Users/fengyanshi15/tmp3/';
dx=0.1;

eta1=load([fdir1 'eta_00001']);

m=length(eta1);
x=(0:m-1)*dx;

show_time=[30:30:160];

clf

[pks,locs]=max(eta1(1,:));
val=sprintf('%.3f',pks*100.);
    text(x(locs),pks,val)
    text(x(locs),pks-pks*0.05,'x 10E-2')

    
  hold on
plot(x,eta1,'r');

grid
xlabel('x (m)')
ylabel('eta (m)')

for num=1:length(show_time)
    fnum=sprintf('%.5d',show_time(num));
    eta1=load([fdir1 'eta_' fnum]);
    [pks,locs]=max(eta1(1,:));
    plot(x,eta1,'b');
    val=sprintf('%.3f',pks*100.);
    text(x(locs),pks,val)
end

print -djpeg dx_p05_ht_p1.jpg


