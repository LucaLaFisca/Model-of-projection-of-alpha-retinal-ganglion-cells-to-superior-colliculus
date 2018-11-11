function [cons,ceq] = optimConstr(x)
b = getGlobalDataTemp;
x_temporal = b(1,:);

c = getGlobalDataSpatial;
x_spatial = c(1,:);

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

%% Spatial filter parameters
dx = x_spatial(2) - x_spatial(1);
x_min = -x_spatial(end); %x-value roughly giving -x border of r.f.
x_max = x_spatial(end);
x_vect = x_min:dx:x_max;

%% Define the spatial kernel
D_x_cen = A_c*exp(-(x_vect.^2)/(2*(sigma_c^2)))/...
    (2*pi*sigma_c^2);
D_x_sur = A_s*exp(-(x_vect.^2)/(2*(sigma_s^2)))/...
    (2*pi*sigma_s^2);

%% temporal filter parameters
tau_vect = x_temporal;

%% Define the temporal kernel
D_t_cen = alpha_cen^2*tau_vect.*exp(-alpha_cen*tau_vect) - ...
   B_cen*beta_cen^2*tau_vect.*exp(-beta_cen*tau_vect);
D_t_sur = alpha_sur^2*tau_vect.*exp(-alpha_sur*tau_vect) - ...
   B_sur*beta_sur^2*tau_vect.*exp(-beta_sur*tau_vect);

%% Constraints
cons = zeros(1,3);
% temporal kernel: positive peak before negative peak
[~,idx_cen_max_temp] = max(D_t_cen);
[~,idx_cen_min_temp] = min(D_t_cen);
cons(1) = idx_cen_max_temp - idx_cen_min_temp;
[~,idx_sur_max_temp] = max(D_t_sur);
[~,idx_sur_min_temp] = min(D_t_sur);
cons(2) = idx_sur_max_temp - idx_sur_min_temp;
% spatial kernel: ensure the right shape for ON/OFF center
[~,idx_max_space] = max(D_x_cen - D_x_sur);
[~,idx_min_space] = min(D_x_cen - D_x_sur);
cons(3) = x_vect(idx_max_space)^2 - x_vect(idx_min_space)^2;

% constraints = cons

cons = max(cons);
ceq = [];