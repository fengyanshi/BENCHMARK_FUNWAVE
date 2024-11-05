clear all

m=3072;
n=1024;
dx=0.5;
x=[0:m-1]*dx;
y=(0:n-1)*dx;
dep=zeros([n,m]);
dep(:,:)=10.0;

dep_shallow=10.0;
dep_deep=10.0;

for j=1:n
    for i=1:m
     dep(j,i)=dep_deep-(dep_deep-dep_shallow)*(j-1.)/(n-1.0);
    end
end

pcolor(x,y,dep),shading flat
colorbar
save -ASCII depth.txt dep