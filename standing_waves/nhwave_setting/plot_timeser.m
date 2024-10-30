fdir1='/Users/fengyanshi/OUTSIDE_Google_1/NHWAVE_IBM/standing_wave/results/layers_6_alpha_0/';
fdir2='/Users/fengyanshi/OUTSIDE_Google_1/NHWAVE_IBM/standing_wave/results/layers_6_alpha_05/';

case1_sta1=load([fdir1 'probe_0001']);
case1_sta2=load([fdir1 'probe_0002']);
case2_sta1=load([fdir2 'probe_0001']);
case2_sta2=load([fdir2 'probe_0002']);

time1=case1_sta1(:,1);
eta1=case1_sta1(:,2);

time2=case2_sta1(:,1);
eta2=case2_sta1(:,2);

figure(1)
clf
plot(time1,eta1,'k',time2,eta2,'r-')
