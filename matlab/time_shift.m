function [ data_shift ] = time_shift( data, mu, sigma, sens, label)
% time shifts every channel of data. Most of the code here is dedicated to
% time shifting each gradiometer pair the same amount
%   Detailed explanation goes here

nChan = size(data,1);
data_shift = zeros(size(data));

%generate distribution
r = round(normrnd(mu,sigma,1,nChan));

% magnetometers
if strcmpi(sens, 'mag')
    for i = 1:nChan
        data_shift(i,:) = circshift(data(i,:), r(i));
    end
else
    % gradiometers - indexed with a 2 or a 3 at the end
    cnt = 0;
    labels = {};
    for i = 1:numel(label)
        curr = label{i};
        if strcmp(curr(end), '2')
            cnt = cnt + 1;
            labels{1,cnt} = label{i};
            labels{2,cnt} = [label{i}(1:end-1), '3'];
        end
    end
    for i = 1:size(labels,2)
        idx2 = strcmp(label, labels{1,i});
        idx3 = strcmp(label, labels{2,i});
        data_shift(idx2,:) = circshift(data(idx2,:), r(i));
        data_shift(idx3,:) = circshift(data(idx3,:), r(i));
    end
end
end

