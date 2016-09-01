clear all;  clc;  close all;

%%%%%%%%%%%%
x0=0.;       x1=53.98;
y0=0.;       y1=0.1;
dx=0.02;   dy=0.02;
%%%%%%%%%%%%
xx = x0 : dx : x1;
yy = y0 : dy : y1;
sizex=length(xx);
sizey=length(yy);
clear x0 x1 y0 y1

[ xgrid, ygrid ] = meshgrid( xx, fliplr(yy) );
dep=zeros( sizey, sizex ) + 0.4;

ind = find( 26 <= xgrid(:) & xgrid(:) <= 32 );
dep(ind) = -(1./20.)*( xgrid(ind)-26. ) + 0.4;
clear ind

ind = find( 32 <= xgrid(:) & xgrid(:) <= 34 );
dep(ind) = 0.1;
clear ind

ind = find( 34 <= xgrid(:) & xgrid(:) <= 37 );
dep(ind) = (1./10.)*( xgrid(ind)-34. ) + 0.1;
clear ind

set( figure(1), 'position', [50 50 800 180] )
surf( xgrid, ygrid, -dep );
shading interp;
view(0,90);  axis tight;
colorbar;

fname = ['dep2D.eps'];
set(gcf, 'PaperPositionMode', 'auto');
print ('-depsc2', fname);
clear fname %xgrid ygrid;

hold off;  

fname=['dep_ny',num2str(size(dep,1)),'mx',num2str(size(dep,2)),'.dat' ]
[FileID]=fopen(  fname ,'wt');
for line=size( dep,1 ):-1:1
    fprintf( FileID, '%10.4f', dep(line,:) );
    fprintf( FileID,'\n');
end; clear FileID
fclose all;