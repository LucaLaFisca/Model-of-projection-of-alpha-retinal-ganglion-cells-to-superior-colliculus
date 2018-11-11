function f = weightsFunc(x)
global temp_step
dt = temp_step;
global reference
final_rates = reference;
global RGCs_rates
save_rates = RGCs_rates;

gain = x(1);
t_shift = x(2);
thresh = x(3);
shift_step = round(t_shift/dt);
weights_size = size(save_rates);
weights_size = weights_size(3)/4;
weights_val = reshape(x(4:end),4,length(x(4:end))/4);
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

shifted_rates = circshift(save_rates,shift_step);
shifted_rates(1:shift_step+1,:,:) = repmat(save_rates(1,:,:),shift_step+1,1,1);
new_rates = reshape(shifted_rates,size(save_rates,1)*size(save_rates,2), size(save_rates,3)/4, 4);
new_weights = reshape(final_weights,size(final_weights,1)*size(final_weights,2),size(final_weights,3));
weighted_rates = zeros(size(new_rates,1),1);
for i = 1:4
    weighted_rates = weighted_rates + new_rates(:,:,i) * new_weights(:,i);
end
computed_rates = reshape(-thresh + gain*weighted_rates,size(save_rates,1),size(save_rates,2));
computed_rates = max(computed_rates,0);


% new_rates = reshape(save_rates,size(save_rates,1)*size(save_rates,2), size(save_rates,3)/4, 4);
% new_weights = reshape(final_weights,size(final_weights,1)*size(final_weights,2),size(final_weights,3));
% computed_rates = zeros(size(new_rates,1),1);
% for i = 1:4
%     computed_rates = computed_rates + new_rates(:,:,i) * new_weights(:,i);
% end
% computed_rates = new_save_rates * new_weights(:);
% thresholded_rates = max(0,computed_rates);
% 
% if isempty(find((computed_rates-thresholded_rates).^2 > 0, 1))
%     disp('no difference after thresholding')
% end

% computed_rates = reshape(computed_rates,size(save_rates,1),size(save_rates,2));

global costFunc
% [val_ref,idx_ref] = max(final_rates);
% [val_comput,idx_comput] = max(computed_rates);
% cmp_max = sum((idx_ref - idx_comput).^4 + (val_ref - val_comput).^2);
% f = sum(sum((final_rates - computed_rates).^2))/size(final_rates,1) + cmp_max;

f = sum(sum((final_rates - computed_rates).^2));
costFunc = f;

%%
% weights_size = size(save_rates);
% weights = zeros(ceil(weights_size(1)/2));
% weights_val = x;
% k=0;
% for j = 1:ceil(weights_size/2)
%     for i = 1:j
%         k = k+1;
%         weights(i,j) = weights_val(k);
%     end
% end
% k=0;
% for i = 1:ceil(weights_size/2)
%     for j = 1:i
%         k = k+1;
%         weights(i,j) = weights_val(k);
%     end
% end
% 
% top_left = weights;
% top_right = fliplr(weights);
% top_right(:,1) = [];
% top = [top_left top_right];
% bottom = flipud(top);
% bottom(1,:) = [];
% 
% final_weights = [top; bottom];


% new_save_rates = reshape(save_rates,size(save_rates,1)*size(save_rates,2), length(t_vect));
% computed_rates = new_save_rates' * final_weights(:);
% for t = 1:length(t_vect)
%     computed_rates(t) = sum(sum(final_weights .* save_rates(:,:,t)));
% end