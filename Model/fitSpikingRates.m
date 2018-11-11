function [B_cen, B_sur, alpha_cen, alpha_sur, beta_cen, beta_sur, A_c, A_s, sigma_c, sigma_s] = ...
    fitSpikingRates(dx, dt, x0,spont_rate,fig_x,fig_t)
r0 = spont_rate;

[x_spatial, y_spatial] = get_spatial_coordinates(dx,fig_x);
[x_temporal, y_temporal] = get_temporal_coordinates(dt,fig_t);

case_ON = strfind(fig_x,'ON_');
if case_ON
    case_ON = true;
    My_spatial_case = 'ON';
    kernel_sign = 1;
else
    case_ON = false;
    My_spatial_case = 'OFF';
    kernel_sign = -1;
end

case_trans = strfind(fig_x,'_trans');
if case_trans
    My_temporal_case = 'transient';
else
    My_temporal_case = 'sustainable';
end 

%% OPTIMIZATION
setGlobalDataSpatial([x_spatial; y_spatial])
setGlobalDataTemp([x_temporal; y_temporal])
setGlobalFunc([r0 case_ON])
fun = @optimFunc;
lb = [.5,.5,1E-2,1E-2,1E-3,1E-3,1E-3,1E-3,1E-3,1E-3];
ub = [1,1,1,1,1,1,1E6,1E6,1E3,1E3];
A = [0,0,-1,0,1,0,0,0,0,0; 0,0,0,-1,0,1,0,0,0,0];
b = [0;0];
Aeq = [];
beq = [];
nonlcon = @optimConstr;
options = optimoptions('fmincon','TolFun',1e-6,'MaxIter',500,'MaxFunEvals',5000);
% options.Display = 'iter';
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
% x = x0;
B_cen = x(1);
B_sur = x(2);
alpha_cen = x(3);
alpha_sur = x(4);
beta_cen = x(5);
beta_sur = x(6);
A_c = x(7);
A_s = x(8);
sigma_c = x(9);
sigma_s = x(10);

%% TEST THE MODEL


%% Spatial filter parameters
dx = x_spatial(2) - x_spatial(1);
x_min = -4*x_spatial(end); %x-value roughly giving -x border of r.f.
x_max = 4*x_spatial(end);
x_vect = x_min:dx:x_max;
[~,center_idx] = max(y_spatial);
surround_idx = length(x_spatial);

%% Define spatial kernel
x_array = repmat(x_vect,size(x_vect,2),1);
D_x_cen = A_c*exp(-(x_array.^2 + (x_array').^2)/(2*(sigma_c^2)))/...
    (2*pi*sigma_c^2);
D_x_sur = A_s*exp(-(x_array.^2 + (x_array').^2)/(2*(sigma_s^2)))/...
    (2*pi*sigma_s^2);

%% temporal filter parameters
dt = x_temporal(2) - x_temporal(1);
tau_vect = x_temporal;
t_vect = x_temporal;

%% Define temporal kernel
D_t_cen = alpha_cen^2*tau_vect.*exp(-alpha_cen*tau_vect) - ...
   B_cen*beta_cen^2*tau_vect.*exp(-beta_cen*tau_vect);
D_t_sur = alpha_sur^2*tau_vect.*exp(-alpha_sur*tau_vect) - ...
   B_sur*beta_sur^2*tau_vect.*exp(-beta_sur*tau_vect);

%% Compute temporal rates for every spot size
spatial_range = 0:dx:40;
[xx, yy] = meshgrid(x_vect,x_vect);
test_stim = zeros(size(xx));
final_integral = zeros(length(spatial_range),length(t_vect)); %set up vector to hold temporal filter values
spatial_rates = zeros(1,length(spatial_range));

i = 0;
for j = spatial_range
    i = i + 1;
    test_stim(((xx).^2+(yy).^2)<=(j/2).^2) = kernel_sign;
    spatial_center_integral = sum(sum(dx^2*D_x_cen.*test_stim));
    spatial_surround_integral = sum(sum(dx^2*D_x_sur.*test_stim));
    total_integral = dt*(kernel_sign * ...
            (D_t_cen*spatial_center_integral - ...
            D_t_sur*spatial_surround_integral));     
    final_integral(i,:) = cumsum(total_integral);
    if length(total_integral) < length(x_temporal)
        final_integral(i,length(x_temporal)-length(total_integral):length(x_temporal)) = sum(total_integral);
    end
    
    peaks = findpeaks(final_integral(i,:));
    if isempty(peaks)
        spatial_rates(i) = r0;
    else
        spatial_rates(i) = r0 + peaks(1);
    end
end
final_rates = r0 + final_integral;
final_rates = max(0,final_rates);
temporal_rates = final_rates(center_idx,:);
temporal_sur_rates = final_rates(surround_idx,:);

%% Plot the rates
My_header = [My_spatial_case ' center ' My_temporal_case];
%temporal
figure
subplot(1,2,1)
plot(x_temporal,y_temporal)
My_title = [My_header ' temporal reference rates'];
title(My_title)
xlabel('time (ms)')
ylabel('spike rate (Hz)')
subplot(1,2,2)
plot(x_temporal,temporal_rates)
My_title = [My_header ' temporal computed rates'];
title(My_title);
xlabel('time (ms)')
ylabel('spike rate (Hz)')
%temporal surround
figure
plot(x_temporal,temporal_sur_rates)
My_title = [My_header ' temporal surround computed rates'];
title(My_title);
xlabel('time (ms)')
ylabel('spike rate (Hz)')
%spatial reference range
figure
subplot(1,2,1)
plot(x_spatial,y_spatial)
My_title = [My_header ' spatial reference rates'];
title(My_title);
xlabel('spot size (degrees)')
ylabel('spike rate (Hz)')
subplot(1,2,2)
plot(x_spatial,spatial_rates(1:length(x_spatial)))
My_title = [My_header ' spatial computed rates'];
title(My_title);
xlabel('spot size (degrees)')
ylabel('spike rate (Hz)')
%spatial complete range
figure
plot(spatial_range, spatial_rates)
My_title = [My_header ' spatial computed rates (big range)'];
title(My_title);
xlabel('spot size (degrees)')
ylabel('spike rate (Hz)')