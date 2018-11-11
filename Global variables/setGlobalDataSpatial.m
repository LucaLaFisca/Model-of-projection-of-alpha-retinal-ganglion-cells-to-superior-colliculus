function setGlobalDataSpatial(val)
global data_spatial
data_spatial = zeros(2,size(val,2));
data_spatial(1,:) = val(1,:);
data_spatial(2,:) = val(2,:);