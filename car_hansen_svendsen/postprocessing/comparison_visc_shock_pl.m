clear all
close all

data=load('ht_setup_031041.txt');
x_data=data(:,1);
ht_data=data(:,2);
set_data=data(:,3);

case1=load('eddy_vis_case.txt');
x1=case1(:,1);
ht1=case1(:,2);
set1=case1(:,3);

case2=load('shock_cap_case.txt');
x2=case2(:,1);
ht2=case2(:,2);
set2=case2(:,3);

subplot(211)
plot(x_data,ht_data,'ko',x1,ht1,'bd',x2,ht2,'rx')
grid
ylabel('wave height (m)')
legend('data','viscosity breaking','shock-capturing breaking','Location','NorthWest')

subplot(212)
plot(x_data,set_data,'ko',x1,set1,'bd',x2,set2,'rx')
grid
ylabel('setup/setdown (m)')
xlabel('distance from toe of beach (m)')