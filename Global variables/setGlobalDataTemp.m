function setGlobalDataTemp(val)
global data_temp
data_temp = zeros(2,size(val,2));
data_temp(1,:) = val(1,:);
data_temp(2,:) = val(2,:);