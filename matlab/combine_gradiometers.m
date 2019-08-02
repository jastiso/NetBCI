function [ data_combined ] = combine_gradiometers( data, ft_flag )
% Get norm of gradiometers data. 

% ft_flag indicates whether the data is formatted as fieldtrip powspctrm,
% or as a psd, output from matlab pwelch method. Those are the only two
% options right now

if ft_flag
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
else
    nGrad = size(data,2); % second dimension of psd gives channels, first gives frequency
    data_combined = zeros(size(data,1),nGrad/2);
    cnt = 1;
    for i = 1:2:nGrad
        data_combined(:,cnt) = sqrt(data(:,i).^2 + data(:,i).^2);
        cnt = cnt + 1;
    end

end

