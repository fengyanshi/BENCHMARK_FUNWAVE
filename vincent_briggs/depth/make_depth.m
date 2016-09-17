clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1=0.;        x2=37.95;
y1=0.;        y2=25;
dx=0.05;      dy=0.05;

xcenter=21.1;   ycenter=12.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = x1 : dx : x2;
y = y1 : dy : y2;
yy=fliplr(y);
[ xgrid, ygrid ]=meshgrid(x',yy'); 
clear x y yy mx ny;

dep=zeros( size(xgrid) )+0.4572;
xr=xgrid-xcenter;  yr=ygrid-ycenter;
ids=find( (xr/3.05).*(xr/3.05) + (yr/3.96).*(yr/3.96) <= 1 );
if( ~isempty(ids) )
    dep(ids) = 0.9144 - 0.762*sqrt( 1 - (xr(ids)/3.81).*(xr(ids)/3.81) - (yr(ids)/4.95).*(yr(ids)/4.95) );
end
clear ids

surf(xgrid, ygrid, dep)
view(0,90)
axis equal; axis tight;
colorbar;
shading interp; colormap jet   

pname = ['vincent_dep.png'];
print ('-dpng', pname);
hold off;  clear pname

fname=['dep_ny',num2str(size(dep,1)),'mx',num2str(size(dep,2)),'.dat' ]
[FileID]=fopen(  fname ,'wt');
for line=size(dep,1):-1:1

    fprintf( FileID, '%12.7f', dep(line,:) );
    fprintf( FileID,'\n');    

end; clear FileID
fclose all;

min(dep(:))
max(dep(:))