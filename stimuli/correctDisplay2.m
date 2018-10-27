function res = correctDisplay2(input, varargin)

if ~isempty(varargin) > 0
    binary = varargin{1};
else
    binary = 0;
end
global Set

% Monitor size and position variables
w = Set.ScreenWidth;  % width of screen, in cm
h = Set.ScreenHeight;  % height of screen, in cm
cx = Set.mouseposx;   % eye x location, in cm
cy = Set.mouseposy; % eye y location, in cm

% Distance to bottom of screen, along the horizontal eye line
zdistBottom = Set.zdistBottom;     % in cm
%zdistTop    = 21;     % in cm

%Shortest distance from screen to eye (defines sf of unmorphed image)




% Alternatively, you can specify the angle of the screen
screenAngle = Set.ScreenAngle;   % in degrees, measured from table surface in front of screen to plane of screen
zdistTop = zdistBottom - (h*sin(deg2rad(90-screenAngle)));

pxXmax = Set.Screenx; % number of pixels in an image that fills the whole screen, x
pxYmax = Set.Screeny; % number of pixels in an image that fills the whole screen, y

% Internal conversions
top = h-cy;
bottom = -cy;
right = cx;
left = cx - w;

% Convert Cartesian to spherical coord
% In image space, x and y are width and height of monitor and z is the
% distance from the eye. I want Theta to correspond to azimuth and Phi to
% correspond to elevation, but these are measured from the x-axis and x-y
% plane, respectively. So I need to exchange the axes this way, prior to
% converting to spherical coordinates:
% orig (image) -> for conversion to spherical coords
% Z -> X
% X -> Y
% Y -> Z
[xi,yi] = meshgrid(1:pxXmax,1:pxYmax);
cart_pointsX = left + (w/pxXmax).*xi;
cart_pointsY = top - (h/pxYmax).*yi;
cart_pointsZ = zdistTop + ((zdistBottom-zdistTop)/pxYmax).*yi;
[~, sphr_pointsPh, ~] ...
            = cart2sph(cart_pointsZ,cart_pointsX,cart_pointsY);
        
        
[~, sphr_pointsPh2, ~] ...
            = cart2sph(cart_pointsZ,cart_pointsY,cart_pointsX); 
%% try a distortion

if binary
   mi = min(min(input));
   ma = max(max(input));
end

%figure out how large input image is and how many deg it spans

[hs, ws ] = size(input);


ranx = ws*Set.RadPerPix;
rany = hs*Set.RadPerPix;

anglesx = -(ranx/2-Set.RadPerPix):Set.RadPerPix:ranx/2;
anglesy = -(rany/2-Set.RadPerPix):Set.RadPerPix:rany/2;

[xang, yang] = meshgrid(anglesx, anglesy*-1);
%r = repmat(10, size(xang, 1), size(xang, 1));
%[z, x, y] = sph2cart(xang, yang, r);

res = interp2(xang, yang, double(input), sphr_pointsPh2, sphr_pointsPh);

if binary
    res = res - min(min(res));
    res = res/max(max(res));
    res = round(res);
    res(res == 1) = ma;
    res(res == 0) = mi;
end

