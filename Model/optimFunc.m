function f = optimFunc(x)
a = getGlobalFunc;
r0 = a(1);
case_ON = a(2);

b = getGlobalDataTemp;
x_temporal = b(1,:);
y_temporal = b(2,:);

c = getGlobalDataSpatial;
x_spatial = c(1,:);
y_spatial = c(2,:);

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

%% ON-center or OFF-center
if case_ON
    kernel_sign = 1;
else
    kernel_sign = -1;
end

%% Spatial filter parameters
dx = x_spatial(2) - x_spatial(1);
x_min = -x_spatial(end); %x-value roughly giving -x border of r.f.
x_max = x_spatial(end);
x_vect = x_min:dx:x_max;

%% Define the spatial kernel
x_array = repmat(x_vect,size(x_vect,2),1);
D_x_cen = A_c*exp(-(x_array.^2 + (x_array').^2)/(2*(sigma_c^2)))/...
    (2*pi*sigma_c^2);
D_x_sur = A_s*exp(-(x_array.^2 + (x_array').^2)/(2*(sigma_s^2)))/...
    (2*pi*sigma_s^2);

%% temporal filter parameters
dt = x_temporal(2) - x_temporal(1);
tau_vect = x_temporal;
t_vect = x_temporal;

%% Define the temporal kernel
D_t_cen = alpha_cen^2*tau_vect.*exp(-alpha_cen*tau_vect) - ...
   B_cen*beta_cen^2*tau_vect.*exp(-beta_cen*tau_vect);
D_t_sur = alpha_sur^2*tau_vect.*exp(-alpha_sur*tau_vect) - ...
   B_sur*beta_sur^2*tau_vect.*exp(-beta_sur*tau_vect);


%% Compute temporal rates for every spot size
[xx, yy] = meshgrid(x_vect,x_vect);
test_stim = zeros(size(xx));
final_integral = zeros(length(x_spatial),length(t_vect)); %set up vector to hold temporal filter values
spatial_rates = zeros(1,length(x_spatial));

i = 0;
for j = x_spatial
    i = i + 1;
    test_stim(((xx).^2+(yy).^2)<=(x_spatial(i)/2).^2) = kernel_sign;
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
        [val_max,idx_max] = max(final_integral(i,:));
        if and(or(idx_max == 1, idx_max == length(final_integral)),...
                ~(val_max == 0))
            spatial_rates(i) = -1E6;
%             disp('error inside empty')
        else
            spatial_rates(i) = r0;
        end
    elseif ~(max(final_integral(i,:)) == peaks(1))
        spatial_rates(i) = -1E6;
%         disp('error outside empty')
    else
        spatial_rates(i) = r0 + peaks(1);
    end
end

final_rates = r0 + final_integral;
final_rates = max(0,final_rates);
[~,center_idx] = max(y_spatial);
temporal_rates = final_rates(center_idx,:);   
% save('RatesInFunc', 'final_rates', 'temporal_rates', 'alpha_cen', 'alpha_sur', 'beta_cen', 'beta_sur', 'A_c', 'A_s', 'sigma_c', 'sigma_s')

%% Cost function temporal
f_temporal = sum((temporal_rates - y_temporal).^2)/length(y_temporal);

%% Cost function positive peaks temporal
[~,idx_max] = max(y_temporal);
[~,idx_computed] = max(temporal_rates);
if abs(idx_max - idx_computed) >= 5
    f_pos_peak_temporal = 1E10;
%     disp('Error peak temporal')
else
    f_pos_peak_temporal = (x_temporal(idx_max) - x_temporal(idx_computed))^2 + ((temporal_rates(idx_computed) - ...
        y_temporal(idx_max)).^2);
end
%% Cost function count positive peaks temporal
[~,idx_space] = max(y_spatial);
[~,idx_max] = max(y_temporal);
[~, idx_flat] = min(abs(y_temporal(idx_max:end) - (min(y_temporal(idx_max:end))+2)));
idx_flat = idx_max-1 + idx_flat;    %because index 0 was idx_max-1 previously
pos_peaks = findpeaks(final_integral(idx_space,1:idx_flat));
f_count_pos_peaks_temporal = length(pos_peaks);

%% Cost function negative count temporal surround
flat = y_temporal(end)-r0;
neg_count = final_integral(final_integral(end,:) < flat/2);
f_neg_count = length(neg_count);

%% Cost function flat temporal
[~, idx_flat] = min(abs(y_temporal(idx_max:end) - (min(y_temporal(idx_max:end))+2)));
idx_flat = idx_max-1 + idx_flat;    %because index 0 was idx_max-1 previously
f_flat_temporal = sum((temporal_rates(idx_flat:end) - ...
    y_temporal(idx_flat:end)).^2)/length(y_temporal(idx_flat:end));

%% Cost function spatial
f_spatial = sum((spatial_rates - y_spatial).^2)/length(y_spatial);

%% Cost function flat spatial
[~,idx_max] = max(y_spatial);
[~, idx_flat] = min(abs(y_spatial(idx_max:end) - (min(y_spatial(idx_max:end))+1)));
idx_flat = idx_max-1 + idx_flat;    %because index 0 was idx_max-1 previously
f_flat_spatial = sum((spatial_rates(idx_flat:end) - y_spatial(idx_flat:end)).^2)/length(y_spatial(idx_flat:end));

%% Cost function peak spatial
[~,idx_max] = max(y_spatial);
[~,idx_computed] = max(spatial_rates);
f_peak_spatial = ((x_spatial(idx_max) - x_spatial(idx_computed))*10)^2 +...
    ((spatial_rates(idx_computed) - y_spatial(idx_max)).^2);

%% Cost function Final
f = f_temporal + f_spatial/2 + f_flat_spatial + f_peak_spatial + f_pos_peak_temporal +...
    f_flat_temporal + f_count_pos_peaks_temporal^1000 + f_neg_count;