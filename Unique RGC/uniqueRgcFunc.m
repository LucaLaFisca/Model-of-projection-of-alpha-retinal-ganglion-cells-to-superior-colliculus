function f = uniqueRgcFunc(x)
global temp_step
dt = temp_step;
global reference
final_rates = reference;
global RGCs_rates
save_rates = RGCs_rates;

gain = x(1);
t_shift = x(2);
spont = x(3);
shift_step = round(t_shift/dt);
weights_val = x(4:end);

final_weights = weights_val;


shifted_rates = circshift(save_rates,[shift_step 0 0]);
shifted_rates(1:shift_step+1,:,:) = repmat(save_rates(1,:,:),shift_step+1,1,1);
new_rates = reshape(shifted_rates,size(save_rates,1)*size(save_rates,2), 4);
new_weights = final_weights;
weighted_rates = zeros(size(new_rates,1),1);
for i = 1:4
    weighted_rates = weighted_rates + new_rates(:,i) * new_weights(i);
end
computed_rates = reshape(spont + gain*weighted_rates,size(save_rates,1),size(save_rates,2));


global costFunc

f = sum(sum((final_rates - computed_rates).^2));
costFunc = f;