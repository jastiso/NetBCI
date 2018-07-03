function [ data_combined ] = combine_gradiometers( data )
% Get norm of gradiometers data from fieldtrip

nGrad = numel(data.label);
data_combined = data;
data_combined.label = [];
data_combined.powspctrm = [];

cnt = 1;
for i = 1:2:nGrad
    data_combined.label{cnt,1} = [data.label{i}, '+', data.label{i+1}(end-3:end)];
    data_combined.powspctrm(:,cnt,:) = sqrt(data.powspctrm(:,i,:).^2 + data.powspctrm(:,i,:).^2);
    cnt = cnt + 1;
end

end

