[x,y] = trace_curve('alpha spikes temporal curve.PNG');
close all
L = 70.7/1000;
xx=0:1/1000:L;
yy=spline(x,y,xx);
figure
plot(xx,yy)

%%
P = 134;
F = 10*(41.5/100);
T = 45.5*(1+41.5/100)/1000;
t=L:1/1000:1;
a = (P-F)*exp(-(t-L)/T)+F;
figure
plot(t,a)

%%
final_x = horzcat(xx, t);
final_y = horzcat(yy, a);
figure
% uiopen('ON_trans_temporal.fig',1)
% hold on
plot(final_x, final_y)
savefig('ON_trans_temporal_surround')