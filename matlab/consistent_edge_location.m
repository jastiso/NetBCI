%% Get lobes of most consistent edges

% Global Variables

addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
addpath(genpath('/Users/stiso/Documents/MATLAB/BCT/'))

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/param/';
save_dir = [top_dir, 'GroupAvg/wpli/analysis/'];
img_dir = [top_dir, 'GroupAvg/wpli/images/'];

% have motor at the end so they will overwrite
regions = [{'Vertex'},{'Left_frontal'}, {'Left_occipital'}, {'Left_parietal'}, {'Left_temporal'}, ...
    {'Right_frontal'}, {'Right_occipital'}, {'Right_parietal'}, {'Right_temporal'}, {'Left_motor'}, {'Right_motor'},...
    {'Frontoparietal'}, {'Parietoccipital'}, {'Frontal'}];

subjs = [1:20];
nSubj = numel(subjs);
nNode = 102;
nEdge = (nNode^2-nNode)/2;
freqs = [7,14;15,30;31,45;55,70];
bands = [{'alpha'},{'beta'},{'low_gamma'}];
sensors = [{'grad'}];

load([save_dir, 'most_consistent_pr.mat'])
high_pr = c_high_max;
high2_pr = c_high2_max;
high3_pr = c_high3_max;
low_pr = c_low_max;

load([save_dir, 'most_consistent.mat'])
high = c_high_max(:,3);
high2 = c_high2_max(:,3);
high3 = c_high3_max(:,3);
low = c_low_max(:,3);

% intialize
high_edges = cell(size(high));
high2_edges = cell(size(high2));
high3_edges = cell(size(high3));
low_edges = cell(size(low));

high_edges_pr = cell(size(high_pr));
high2_edges_pr = cell(size(high2_pr));
high3_edges_pr = cell(size(high3_pr));
low_edges_pr = cell(size(low_pr));

for j = 1:numel(bands)
    % initialize
    high_curr = high{j};
    high2_curr = high2{j};
    high3_curr = high3{j};
    low_curr = low{j};
    
    high_edges{j} = repmat({''},size(high_curr));
    high2_edges{j} = repmat({''},size(high2_curr));
    high3_edges{j} = repmat({''},size(high3_curr));
    low_edges{j} = repmat({''},size(low_curr));
    
    high_curr_pr = high_pr{j};
    high2_curr_pr = high2_pr{j};
    high3_curr_pr = high3_pr{j};
    low_curr_pr = low_pr{j};
    
    high_edges_pr{j} = repmat({''},size(high_curr_pr));
    high2_edges_pr{j} = repmat({''},size(high2_curr_pr));
    high3_edges_pr{j} = repmat({''},size(high3_curr_pr));
    low_edges_pr{j} = repmat({''},size(low_curr_pr));
    
    for i = 1:numel(regions)
        load([top_dir, 'montages/', regions{i}, '_idx.mat'])
        
        % empirical
        % loop through edges in high
        for m = 1:size(high_curr,1)
            if idx(high_curr(m,1))
                high_edges{j}{m,1} = regions{i};
            end
            if idx(high_curr(m,2))
                high_edges{j}{m,2} = regions{i};
            end
        end
        high_edges{j}(:,3) = num2cell(high_curr(:,3));
        
        % loop through edges in high2
        for m = 1:size(high2_curr,1)
            if idx(high2_curr(m,1))
                high2_edges{j}{m,1} = regions{i};
            end
            if idx(high2_curr(m,2))
                high2_edges{j}{m,2} = regions{i};
            end
        end
        high2_edges{j}(:,3) = num2cell(high2_curr(:,3));
        
        % loop through edges in high3
        for m = 1:size(high3_curr,1)
            if idx(high3_curr(m,1))
                high3_edges{j}{m,1} = regions{i};
            end
            if idx(high3_curr(m,2))
                high3_edges{j}{m,2} = regions{i};
            end
        end
        high3_edges{j}(:,3) = num2cell(high3_curr(:,3));
        
        % loop through edges in low
        for m = 1:size(low_curr,1)
            if idx(low_curr(m,1))
                low_edges{j}{m,1} = regions{i};
            end
            if idx(low_curr(m,2))
                low_edges{j}{m,2} = regions{i};
            end
        end
        low_edges{j}(:,3) = num2cell(low_curr(:,3));
        
        
        
        
        % UPR
        % loop through edges in high
        for m = 1:size(high_curr_pr,1)
            if idx(high_curr_pr(m,1))
                high_edges_pr{j}{m,1} = regions{i};
            end
            if idx(high_curr_pr(m,2))
                high_edges_pr{j}{m,2} = regions{i};
            end
        end
        high_edges_pr{j}(:,3) = num2cell(high_curr_pr(:,3));
        
        % loop through edges in high2
        for m = 1:size(high2_curr_pr,1)
            if idx(high2_curr_pr(m,1))
                high2_edges_pr{j}{m,1} = regions{i};
            end
            if idx(high2_curr_pr(m,2))
                high2_edges_pr{j}{m,2} = regions{i};
            end
        end
        high2_edges_pr{j}(:,3) = num2cell(high2_curr_pr(:,3));
        
        % loop through edges in high3
        for m = 1:size(high3_curr_pr,1)
            if idx(high3_curr_pr(m,1))
                high3_edges_pr{j}{m,1} = regions{i};
            end
            if idx(high3_curr_pr(m,2))
                high3_edges_pr{j}{m,2} = regions{i};
            end
        end
        high3_edges_pr{j}(:,3) = num2cell(high3_curr_pr(:,3));
        
        % loop through edges in low
        for m = 1:size(low_curr_pr,1)
            if idx(low_curr_pr(m,1))
                low_edges_pr{j}{m,1} = regions{i};
            end
            if idx(low_curr_pr(m,2))
                low_edges_pr{j}{m,2} = regions{i};
            end
        end
        low_edges_pr{j}(:,3) = num2cell(low_curr_pr(:,3));
        
    end
    
    % sort each row
     high_edges{j}(:,1:2) = sort(high_edges{j}(:,1:2)')';
     high2_edges{j}(:,1:2) = sort(high2_edges{j}(:,1:2)')';
     high3_edges{j}(:,1:2) = sort(high3_edges{j}(:,1:2)')';
     low_edges{j}(:,1:2) = sort(low_edges{j}(:,1:2)')';
     high_edges_pr{j}(:,1:2) = sort(high_edges_pr{j}(:,1:2)')';
     high2_edges_pr{j}(:,1:2) = sort(high2_edges_pr{j}(:,1:2)')';
     high3_edges_pr{j}(:,1:2) = sort(high3_edges_pr{j}(:,1:2)')';
     low_edges_pr{j}(:,1:2) = sort(low_edges_pr{j}(:,1:2)')';
    
     % order rows by occurance
    high_edges{j} = sortrows(high_edges{j},3);
    high2_edges{j} = sortrows(high2_edges{j},3);
    high3_edges{j} = sortrows(high3_edges{j},3);
    low_edges{j} = sortrows(low_edges{j},3);
    high_edges_pr{j} = sortrows(high_edges_pr{j},3);
    high2_edges_pr{j} = sortrows(high2_edges_pr{j},3);
    high3_edges_pr{j} = sortrows(high3_edges_pr{j},3);
    low_edges_pr{j} = sortrows(low_edges_pr{j},3);
    
    % change to a table so I can sort things later
    high_edges{j} = table(high_edges{j}(:,1),high_edges{j}(:,2),high_edges{j}(:,3));
    high2_edges{j} = table(high2_edges{j}(:,1),high2_edges{j}(:,2),high2_edges{j}(:,3));
    high3_edges{j} = table(high3_edges{j}(:,1),high3_edges{j}(:,2),high3_edges{j}(:,3));
    low_edges{j} = table(low_edges{j}(:,1),low_edges{j}(:,2),low_edges{j}(:,3));
    
    high_edges_pr{j} = table(high_edges_pr{j}(:,1),high_edges_pr{j}(:,2),high_edges_pr{j}(:,3));
    high2_edges_pr{j} = table(high2_edges_pr{j}(:,1),high2_edges_pr{j}(:,2),high2_edges_pr{j}(:,3));
    high3_edges_pr{j} = table(high3_edges_pr{j}(:,1),high3_edges_pr{j}(:,2),high3_edges_pr{j}(:,3));
    low_edges_pr{j} = table(low_edges_pr{j}(:,1),low_edges_pr{j}(:,2),low_edges_pr{j}(:,3));
end

% get unique edges
for j = 1:numel(bands)
    [high_edges_unique{j}, a, c] = unique(high_edges{j}(:,1:2), 'rows', 'stable');
    occurences = zeros(size(high_edges_unique{j},1),1);
    for i = 1:numel(occurences)
        occurences(i) = numel(find(c == i));
    end
    high_edges_unique{j} = [high_edges_unique{j}, table(occurences, 'VariableNames',{'occurances'}), high_edges{j}(a,3)];
    
    [high2_edges_unique{j}, a, c] = unique(high2_edges{j}(:,1:2), 'rows', 'stable');
    occurences = zeros(size(high2_edges_unique{j},1),1);
    for i = 1:numel(occurences)
        occurences(i) = numel(find(c == i));
    end
    high2_edges_unique{j} = [high2_edges_unique{j}, table(occurences, 'VariableNames',{'occurances'}), high2_edges{j}(a,3)];
    
    [high3_edges_unique{j}, a, c] = unique(high3_edges{j}(:,1:2), 'rows', 'stable');
    occurences = zeros(size(high3_edges_unique{j},1),1);
    for i = 1:numel(occurences)
        occurences(i) = numel(find(c == i));
    end
    high3_edges_unique{j} = [high3_edges_unique{j}, table(occurences, 'VariableNames',{'occurances'}), high3_edges{j}(a,3)];
    
    [low_edges_unique{j}, a, c] = unique(low_edges{j}(:,1:2), 'rows', 'stable');
    occurences = zeros(size(low_edges_unique{j},1),1);
    for i = 1:numel(occurences)
        occurences(i) = numel(find(c == i));
    end
    low_edges_unique{j} = [low_edges_unique{j}, table(occurences, 'VariableNames',{'occurances'}), low_edges{j}(a,3)];
    
    
    [high_edges_unique_pr{j}, a, c] = unique(high_edges_pr{j}(:,1:2), 'rows', 'stable');
    occurences = zeros(size(high_edges_unique_pr{j},1),1);
    for i = 1:numel(occurences)
        occurences(i) = numel(find(c == i));
    end
    high_edges_unique_pr{j} = [high_edges_unique_pr{j}, table(occurences, 'VariableNames',{'occurances'}), high_edges_pr{j}(a,3)];
    
    [high2_edges_unique_pr{j}, a, c] = unique(high2_edges_pr{j}(:,1:2), 'rows', 'stable');
    occurences = zeros(size(high2_edges_unique_pr{j},1),1);
    for i = 1:numel(occurences)
        occurences(i) = numel(find(c == i));
    end
    high2_edges_unique_pr{j} = [high2_edges_unique_pr{j}, table(occurences, 'VariableNames',{'occurances'}), high2_edges_pr{j}(a,3)];
    
    [high3_edges_unique_pr{j}, a, c] = unique(high3_edges_pr{j}(:,1:2), 'rows', 'stable');
    occurences = zeros(size(high3_edges_unique_pr{j},1),1);
    for i = 1:numel(occurences)
        occurences(i) = numel(find(c == i));
    end
    high3_edges_unique_pr{j} = [high3_edges_unique_pr{j}, table(occurences, 'VariableNames',{'occurances'}), high3_edges_pr{j}(a,3)];
    
    [low_edges_unique_pr{j}, a, c] = unique(low_edges_pr{j}(:,1:2), 'rows', 'stable');
    occurences = zeros(size(low_edges_unique_pr{j},1),1);
    for i = 1:numel(occurences)
        occurences(i) = numel(find(c == i));
    end
    low_edges_unique_pr{j} = [low_edges_unique_pr{j}, table(occurences, 'VariableNames',{'occurances'}), low_edges_pr{j}(a,3)];
end


save([save_dir, 'consistent_edges.mat'], 'high_edges_unique', 'high2_edges_unique', 'high3_edges_unique', 'low_edges_unique', ...
    'high_edges_unique_pr', 'high2_edges_unique_pr', 'high3_edges_unique_pr', 'low_edges_unique_pr')

