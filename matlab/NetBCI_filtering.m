function [ data_filt ] = NetBCI_filtering( data )
%Filters mEG data
% expects fieldtrip format
% filter
% high pass filter
nTrial = size(data.trial{1},2);
for i = 1:nTrial
    [a,b] = butter(4, .1/(data.fsample/2), 'high');
    data.trial{i} = filtfilt(a,b,data.trial{i}')';
end
data_filt = data;

% low pass filter
%                 [a,b] = butter(4, 160/(data.fsample/2), 'low');
%                 data.trial{1} = filtfilt(a,b,data.trial{1}')';


end

