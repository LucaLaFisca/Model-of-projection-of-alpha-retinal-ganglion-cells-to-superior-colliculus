function [gain, t_shift, spont, weights_val] = uniqueRgcWeights(save_rates,final_rates)
% global time_ref
% time_ref = t_vect;
global reference
reference = final_rates;
global RGCs_rates
RGCs_rates = save_rates;
global costFunc

%%
My_struct = struct('costFunc',1, 'gain',1, 't_shift',20, 'spont',40, 'weights',ones(1,4));
% Final_struct = struct('costFunc',1, 'gain',1, 't_shift',20, 'weights',ones(4,nb_weights));

fun = @uniqueRgcFunc;
lb = [0 20 0 zeros(1,4)];
ub = [1E6 80 1E2 ones(1,4)];
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];

disp('In first loop')
options = optimoptions('fmincon','TolFun',1e-6,'MaxIter',5,'MaxFunEvals',1000);
for i = 1:20
    x0 = [rand(1)*1E2, rand(1)*60+20, rand(1)*1E2 rand(1,4)];
    x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    f = costFunc;
    current_gain = x(1);
    current_t_shift = x(2);
    current_spont = x(3);
    current_weights = x(4:end);
    My_struct(i) = struct('costFunc',f, 'gain',current_gain, 't_shift',...
        current_t_shift, 'spont',current_spont, 'weights',current_weights);
end
disp('In second loop')
My_table = struct2table(My_struct);
func_list = My_table.costFunc;
[~,idx] = min(func_list);
gain = My_table.gain(idx);
t_shift = My_table.t_shift(idx);
spont = My_table.spont(idx);
weights_val = My_table.weights(idx,:);
x0 = [gain t_shift spont weights_val];
options = optimoptions('fmincon','TolFun',1e-6,'MaxIter',500,'MaxFunEvals',5000);
% options.Display = 'iter';
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
f = costFunc;
current_gain = x(1);
current_t_shift = x(2);
current_spont = x(3);
current_weights = x(4:end);
Final_struct = struct('costFunc',f, 'gain',current_gain, 't_shift',current_t_shift, 'spont',current_spont, 'weights',current_weights);
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
spont = Final_struct.spont;
weights_val = Final_struct.weights;
% weights = zeros(ceil(sqrt(weights_size)/2),ceil(sqrt(weights_size)/2),4);
% 
% 
% k=0;
% for j = 1:ceil(sqrt(weights_size)/2)
%     for i = 1:j
%         k = k+1;
%         weights(i,j,:) = weights_val(:,k);
%     end
% end
% k=0;
% for i = 1:ceil(sqrt(weights_size)/2)
%     for j = 1:i
%         k = k+1;
%         weights(i,j,:) = weights_val(:,k);
%     end
% end
% 
% top_left = weights;
% top_right = fliplr(weights);
% top_right(:,1,:) = [];
% top = [top_left top_right];
% bottom = flipud(top);
% bottom(1,:,:) = [];

% final_weights = [top; bottom];