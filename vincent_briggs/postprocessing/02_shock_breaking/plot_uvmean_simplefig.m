clear all; clc;

%%%%%%%%%%%%%%%%%%
mx=760;     ny=501;
dx=0.05;     dy=0.05;

ufname='../result_01/umean_01000';
vfname='../result_01/vmean_01000';
depname='../depth/dep_ny501mx760.dat';
%%%%%%%%%%%%%%%%%%

x = 0 : dx : (mx-1)*dx;
y = 0 : dy : (ny-1)*dy;
x=x';

[xgrid,ygrid]=meshgrid(x,y);  clear x y;

um=load( ufname );
vm=load( vfname );
dep=load( depname );

set(figure(1),'position',[50 50 365 550])

sk=10;
quiver( xgrid(1:sk:end, 1:sk:end), ygrid(1:sk:end, 1:sk:end), ...
    um(1:sk:end, 1:sk:end), vm(1:sk:end, 1:sk:end), 'k', 'LineWidth', 1.5, ...
    'MaxHeadSize', 5.0)

hold on
contour( xgrid, ygrid, dep, [0.1 0.4572], 'r--', 'LineWidth', 2.5 )

xlim([15 34])
ylim([0 25])

xlabel( 'X (m)', 'fontsize', 11 );
ylabel( 'Y (m)', 'fontsize', 11 );