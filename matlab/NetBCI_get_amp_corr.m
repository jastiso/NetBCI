function [ amp_corr ] = NetBCI_get_amp_corr( data, t, srate, f_range )
%Get amplitude cenvelope correlations from data in fieldtrip format

nTrial = numel(data.trial);
nNode = size(data.trial{1},1);
%h = zeros(nTrial, nNode, numel(t)-1);
h = [];
for m = 1:nTrial
    %[a,b] = butter(4, [f_range(1)/(srate/2), f_range(end)/(srate/2)] );
    %h(m,:,:) = filtfilt(a,b,chunk_mag.trial{m}')';
    %h(m,:,:) = eegfilt(data.trial{m}, srate, f_range(1), f_range(2), 0, [], [], 'fir1');
    curr = data.trial{m};
    h(m,:,:) = abs(my_hilbert(curr, srate, f_range(1), f_range(2)));
    %[a,b] = butter(4, 1/(srate/2), 'low' );
    %h(m,:,:) = filtfilt(a,b,squeeze(h(m,:,:)));
end
h = squeeze(mean(h));
% get amplitude corr
amp_corr = corr(h');
amp_corr = amp_corr - eye(nNode);


end

