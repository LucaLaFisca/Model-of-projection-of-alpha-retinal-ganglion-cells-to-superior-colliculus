function [gain, t_shift, thresh, final_weights] = computeWeights(save_rates,final_rates)
% global time_ref
% time_ref = t_vect;
global reference
reference = final_rates;
global RGCs_rates
RGCs_rates = save_rates;
global costFunc

%%

weights_size = size(save_rates);
weights_size = weights_size(3)/4;
weights_vect = zeros(1,ceil(sqrt(weights_size)/2));
nb_weights = sum(1:size(weights_vect,2));

My_struct = struct('costFunc',1, 'gain',1, 't_shift',20, 'thresh',40, 'weights',ones(1,4*nb_weights));
% Final_struct = struct('costFunc',1, 'gain',1, 't_shift',20, 'weights',ones(4,nb_weights));

fun = @weightsFunc;
% lb = zeros(4,nb_weights);
% ub = ones(4,nb_weights);
lb = [0 20 0 zeros(1,4*nb_weights)];
ub = [1E6 80 1E2 ones(1,4*nb_weights)];
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];

% for k = 1:2
disp('In first loop')
options = optimoptions('fmincon','TolFun',1e-6,'MaxIter',5,'MaxFunEvals',1000);
%     options.Display = 'iter';
for i = 1:20
    x0 = [rand(1)*1E2, rand(1)*60+20, rand(1)*1E2 rand(1,4*nb_weights)];
    x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    f = costFunc;
    current_gain = x(1);
    current_t_shift = x(2);
    current_thresh = x(3);
%         current_weights = reshape(x(3:end),4,length(x(3:end))/4);
    current_weights = x(4:end);
    My_struct(i) = struct('costFunc',f, 'gain',current_gain, 't_shift',...
        current_t_shift, 'thresh',current_thresh, 'weights',current_weights);
end
disp('In second loop')
My_table = struct2table(My_struct);
func_list = My_table.costFunc;
[~,idx] = min(func_list);
gain = My_table.gain(idx);
t_shift = My_table.t_shift(idx);
thresh = My_table.thresh(idx);
weights_val = My_table.weights(idx,:);
x0 = [gain t_shift thresh weights_val];
options = optimoptions('fmincon','TolFun',1e-6,'MaxIter',500,'MaxFunEvals',5000);
% options.Display = 'iter';
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
f = costFunc;
current_gain = x(1);
current_t_shift = x(2);
current_thresh = x(3);
current_weights = reshape(x(4:end),4,length(x(4:end))/4);
Final_struct = struct('costFunc',f, 'gain',current_gain, 't_shift',current_t_shift, 'thresh',current_thresh, 'weights',current_weights);
% end

%%
% Final_table = struct2table(Final_struct)
% func_list = Final_table.costFunc;
% [~,idx] = min(func_list);
% gain = Final_table.gain(idx)
% t_shift = Final_table.t_shift(idx)
% weights_val = cell2mat(Final_table.weights(idx))

func_list = Final_struct.costFunc
gain = Final_struct.gain;
t_shift = Final_struct.t_shift;
thresh = Final_struct.thresh;
weights_val = Final_struct.weights;
weights = zeros(ceil(sqrt(weights_size)/2),ceil(sqrt(weights_size)/2),4);


k=0;
for j = 1:ceil(sqrt(weights_size)/2)
    for i = 1:j
        k = k+1;
        weights(i,j,:) = weights_val(:,k);
    end
end
k=0;
for i = 1:ceil(sqrt(weights_size)/2)
    for j = 1:i
        k = k+1;
        weights(i,j,:) = weights_val(:,k);
    end
end

top_left = weights;
top_right = fliplr(weights);
top_right(:,1,:) = [];
top = [top_left top_right];
bottom = flipud(top);
bottom(1,:,:) = [];

final_weights = [top; bottom];






%% Load onset/offset time vectors
% % my_path = input('Enter path of "EVENT" data (from "recordings" directory)\n' ,'s');
% % EVENT = load(my_path);
% % stim_onset = EVENT.strons.Stim*1000;
% % stim_offset = EVENT.strons.Targ*1000;
% 
% 
% 
% 
% 
% 
% %%
% t = 0;
% t_vect = px;
% 
% save_rates = zeros(length(t_vect),size(Mua,2),(2*nb_RGCs+1)*(2*nb_RGCs+1));
% save_rates(1:length(t_vect<0),:,:) = spont_rate;
% 
% i=0;
% for k1 = -nb_RGCs:nb_RGCs
%     i = i+1;
%     j = 0;
%     for k2 = -nb_RGCs : nb_RGCs
%         j = j+1;
%         for trial = 1:size(Mat,1)
% 
%             RGC_stim = stim(xx + k1*alpha_size, yy + k2*alpha_size, trial);
%             figure
%             surface(xx,yy,stim(:,:,trial)); colormap gray
%             title('center')
%             figure
%             surface(xx,yy,RGC_stim(:,:,trial)); colormap gray
%             title('test')
%             pause()
%             spatial_center_integral = sum(sum(dx^2*D_x_cen.*RGC_stim));
%             spatial_surround_integral = sum(sum(dx^2*D_x_sur.*RGC_stim));
% 
%             total_integral = dt*(kernel_sign * ...
%                     (D_t_cen*spatial_center_integral - ...
%                     D_t_sur*spatial_surround_integral));
% 
%             final_integral(i,:) = cumsum(total_integral);
%             if length(total_integral) < length(x_temporal)
%                 final_integral(i,length(x_temporal)-length(total_integral):length(x_temporal)) = sum(total_integral);
%             end
%         end
%     end
% end
% 
% for trial = 1:size(Mat,1)
% %     dt = (stim_offset(trial)-stim_onset(trial))/period;
% %     t_vect = stim_onset(trial):dt:stim_offset(trial);
%     i=0;
%     for k1 = -nb_RGCs:nb_RGCs
%         i = i+1;
%         j = 0;
%         for k2 = -nb_RGCs : nb_RGCs
%             j = j+1;
%             Target_t_vect_long = zeros(size(xx,1),size(xx,1),length(tau_vect));
%             stim = zeros(size(xx));
%             for t_trial = 1:length(t_vect)
%                 t = t+1;
%                 Target_t_vect_long = circshift(Target_t_vect_long,[0 0 1]);
%                 k = find(My_stim.time <= t_vect(t), 1, 'last');
%                 if ~(My_stim.sprite(k) == 0)
%                     stim = zeros(size(xx));
%                     stim(((xx-My_stim.centerx(k) + k1*alpha_size).^2 + ...
%                         (yy-My_stim.centery(k) + k2*alpha_size).^2) <= ...
%                         (My_stim.size(k)/2)^2) = My_stim.sprite(k);
%                     %"(x+xc)^2" instead of "(x-xc)^2" because we
%                     %simulate alpha position from stimulus motion
%                     stim(xx > screenx/2-stim_centerx - k1*alpha_size | ...
%                         xx < -screenx/2-stim_centerx - k1*alpha_size | ...
%                         yy > screeny/2-stim_centery - k2*alpha_size | ...
%                         yy < -screeny/2-stim_centery - k2*alpha_size) = 0;
%                     %Limit the stimulus respecting the size of the screen
%                 end
%                 Target_t_vect_long(:,:,1) = stim;
%                 spike_rate = spont_rate + sum(dt*sum(sum(dx^2*D_x_t_ON_trans.*...
%                     Target_t_vect_long)));
%                 %implement the threshold
%                 spike_rate = max(0,spike_rate);
%                 save_rates(i,j,t) = spike_rate;
%             end
%         end
%     end
% end
