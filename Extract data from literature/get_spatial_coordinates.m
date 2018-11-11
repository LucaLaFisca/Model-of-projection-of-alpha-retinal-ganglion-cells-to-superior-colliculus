function [x,y] = get_spatial_coordinates(dx,fig)
%%Get coordinates from figure
open(fig);
h = gcf;
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children');
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');
step = 1;
while(xdata(1+step)-xdata(1) < dx)
    step = step + 1;
end
x = xdata(1:step:end);
y = ydata(1:step:end);
close
