function [x,y] = get_temporal_coordinates(dt,fig)
%%Get coordinates from figure
open(fig);
title('curve to fit')
% open('ON_transient_temporal.fig')
h = gcf;
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children');
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');
step = 1;
while(xdata(1+step)-xdata(1) < dt/1000)
    step = step + 1;
end
% xdata = xdata(xdata < 0.8);
% ydata = ydata(xdata < 0.8);
x = xdata(1:step:end)*1000;     %1000 converts to ms
y = ydata(1:step:end);
close