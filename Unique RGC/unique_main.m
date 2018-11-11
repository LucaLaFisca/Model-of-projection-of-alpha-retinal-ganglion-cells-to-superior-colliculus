%% Can be launch only after the main if the workspace wasn't cleared

unique_RGC_RFmap = zeros(size(mean_RGC_RFmap,1),size(mean_RGC_RFmap,2),4);
wrong_unique_RGC = zeros(size(mean_RGC_RFmap,1),size(mean_RGC_RFmap,2),4);
central_pos = [25 25+49 25+2*49 25+3*49];
unique_RGC_RFmap = mean_RGC_RFmap(:,:,central_pos);
wrong_unique_RGC = flipud(unique_RGC_RFmap);

total_time = tic;
[gain, t_shift, spont, weights] = uniqueRgcWeights(unique_RGC_RFmap,RFmapMUA(:,:,19));
shift_step = round(t_shift/dt);

[wrong_gain, wrong_t_shift, wrong_spont, wrong_weights] = computeWeights(wrong_unique_RGC,RFmapMUA(:,:,19));
wrong_shift_step = round(wrong_t_shift/dt);

toc(total_time)

%% Final rates for RFmap
shifted_rates = circshift(unique_RGC_RFmap,[shift_step 0 0]);
shifted_rates(1:shift_step+1,:,:) = repmat(unique_RGC_RFmap(1,:,:),shift_step+1,1,1);
new_rates = reshape(shifted_rates,size(unique_RGC_RFmap,1)*size(unique_RGC_RFmap,2), 4);
new_weights = weights;
weighted_rates = zeros(size(new_rates,1),1);
for i = 1:4
    weighted_rates = weighted_rates + new_rates(:,i) * new_weights(i);
end
unique_good_rates = reshape(spont + gain*weighted_rates,size(unique_RGC_RFmap,1),size(unique_RGC_RFmap,2));

%% Final rates for wrong RFmap
shifted_rates = circshift(wrong_RFmap_rates,[wrong_shift_step 0 0]);
shifted_rates(1:wrong_shift_step+1,:,:) = repmat(wrong_RFmap_rates(1,:,:),wrong_shift_step+1,1,1);
new_rates = reshape(shifted_rates,size(wrong_RFmap_rates,1)*size(wrong_RFmap_rates,2), size(wrong_RFmap_rates,3)/4, 4);
new_weights = reshape(wrong_weights,size(wrong_weights,1)*size(wrong_weights,2),size(wrong_weights,3));
weighted_rates = zeros(size(new_rates,1),1);
for i = 1:4
    weighted_rates = weighted_rates + new_rates(:,i) * new_weights(i);
end
unique_wrong_rates = reshape(wrong_spont + wrong_gain*weighted_rates,size(wrong_RFmap_rates,1),size(wrong_RFmap_rates,2));
