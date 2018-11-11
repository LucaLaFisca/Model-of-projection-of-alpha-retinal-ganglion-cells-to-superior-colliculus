% close all
%% Test parameters in the book
dx = .7;
x_vect = -30:dx:30;
dt = 1000/60;
tau_vect = 0:dt:1000;
kernel_sign = 1;
spont_rate = 2.5;

param = [.967 .96699 1/26 1/32 1/48 1/60 1200 700 .93 4.3];
B_cen = param(1);
B_sur = param(2);
alpha_cen = param(3);
alpha_sur = param(4);
beta_cen = param(5);
beta_sur = param(6);
A_c = param(7);
A_s = param(8);
sigma_c = param(9);
sigma_s = param(10);

% B_cen = 0.83;
% B_sur = 0.8;
% alpha_cen = 1/17;
% alpha_sur = 1/18;
% beta_cen = 1/45;
% beta_sur = 1/40;
% A_c = 580;
% A_s = 385;
% sigma_c = 1;
% sigma_s = 3;

%% DISPLAY SPATIAL FILTER 2D
D_x_cen = A_c*exp(-(x_vect.^2)/...
    (2*(sigma_c^2)))/(2*pi*sigma_c^2);
D_x_sur = A_s*exp(-(x_vect.^2)/...
    (2*(sigma_s^2)))/(2*pi*sigma_s^2);
D_x = kernel_sign * (D_x_cen - D_x_sur);

% figure
% plot(x_vect, D_x)
% title('spatial filter 2D')

%% DESIGN TEMPORAL FILTER
D_t_cen = alpha_cen^2*tau_vect.*exp(-alpha_cen*tau_vect) - ...
   B_cen*beta_cen^2*tau_vect.*exp(-beta_cen*tau_vect);
D_t_sur = alpha_sur^2*tau_vect.*exp(-alpha_sur*tau_vect) - ...
   B_sur*beta_sur^2*tau_vect.*exp(-beta_sur*tau_vect);

% figure
% subplot(2,1,1)
% plot(tau_vect,D_t_cen)
% title('D-t-cen-trans')
% subplot(2,1,2)
% plot(tau_vect,D_t_sur)
% title('D-t-sur-trans')

%% Test spatio_temp kernel of book

D_x_t = zeros(length(x_vect),length(tau_vect));
for k = 1:length(tau_vect)
    D_x_t(:,k) = kernel_sign * (D_t_cen(k)*D_x_cen - D_t_sur(k)*D_x_sur);
end
% figure
% surf(D_x_t)
% title('spatio-temp filter')

%% DESIGN THE SPATIAL FILTER (kernel D_x)
x_array = repmat(x_vect,size(x_vect,2),1);
D_x_cen = A_c*exp(-(x_array.^2 + (x_array').^2)/...
    (2*(sigma_c^2)))/(2*pi*sigma_c^2);
D_x_sur = A_s*exp(-(x_array.^2 + (x_array').^2)/...
    (2*(sigma_s^2)))/(2*pi*sigma_s^2);
D_x = kernel_sign * (D_x_cen - D_x_sur);
% figure
% surf(D_x)
% title('Spatial kernel 3D')

%% Faster way
diam_vect = 0:.2:30;
r0 = spont_rate;
t_vect = tau_vect;
% dt = tau_vect(2) - tau_vect(1);
% dx = x_vect(end) - x_vect(end-1);
[xx, yy] = meshgrid(x_vect,x_vect);
test_stim = zeros(size(xx));
final_integral = zeros(length(diam_vect),length(t_vect)); %set up vector to hold temporal filter values
spatial_rates = zeros(1,length(diam_vect));
i = 0;
for j = diam_vect
    i = i + 1;
    test_stim(((xx).^2+(yy).^2)<=(j/2).^2) = kernel_sign;

    spatial_center_integral = sum(sum(dx^2*D_x_cen.*test_stim));
    spatial_surround_integral = sum(sum(dx^2*D_x_sur.*test_stim));
    
    total_integral = dt * (kernel_sign * ...
            (D_t_cen*spatial_center_integral - ...
            D_t_sur*spatial_surround_integral));     
    final_integral(i,:) = cumsum(total_integral);
    
    peaks = findpeaks(final_integral(i,:));
    if isempty(peaks)
        spatial_rates(i) = r0;
    else
        spatial_rates(i) = r0 + peaks(1);
    end
end
final_rates = r0 + final_integral;
[~,idx_max_spatial] = max(spatial_rates);

%%Display
figure
plot(diam_vect,spatial_rates)
title('computed spatial rates');

figure
plot(t_vect,final_rates(idx_max_spatial,:))
title('final integral center')

figure
plot(t_vect,final_rates(length(diam_vect),:))
title('final integral surround')
