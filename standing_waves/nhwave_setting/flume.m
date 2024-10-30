clear all
m=200;
n=10;
dep=zeros([n,m]);
dep(:,:)=2;
slope=0.1;
for i=m-25:m
    dep(:,i)=2.0-(i-m+25)*slope;
end

save -ASCII wave_flume.txt dep
