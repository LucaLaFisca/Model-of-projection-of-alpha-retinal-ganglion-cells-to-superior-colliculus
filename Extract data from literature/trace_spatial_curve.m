[x,y] = trace_curve('alpha spikes curves.PNG');

xx=0:800;
% xx = 0:dx:retina_to_angle(1.6);
% x_test = retina_to_angle(x*2);
yy=spline(x,y,xx);
% yy = interp1(x,y,xx,'spline');
xx = retina_to_angle(xx*2/1000);
figure
plot(xx,yy)
savefig('ON_transient')

