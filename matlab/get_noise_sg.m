%% Find noise sungraphs
% for each band, get the indices of the noise sugbraphs

addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
R_dir = '/Users/stiso/Documents/R/NetBCI/data/wpli/';

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/';
save_dir = [top_dir, 'GroupAvg/wpli/analysis/'];
img_dir = [top_dir, 'GroupAvg/wpli/images/'];
% make directories
if ~exist(img_dir, 'dir')
    mkdir(img_dir)
end
% make directories
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end
% make directories
if ~exist(R_dir, 'dir')
    mkdir(R_dir)
end

Subj = [1:20];
nSubj = numel(Subj);
freqs = [7,14;15,30;31,45;55,70];
bands = [{'alpha'}, {'beta'}, {'low_gamma'}, {'gamma'}];
sensors = [{'grad'}];
regions = [{'Left_frontal'}, {'Left_occipital'}, {'Left_parietal'}, {'Left_temporal'}, ...
    {'Right_frontal'}, {'Right_occipital'}, {'Right_parietal'}, {'Right_temporal'}, {'Vertex'}];

nNode = 102;
nEdges = (nNode^2-nNode)/2;

% initialize
noise_sg = cell(numel(bands),nSubj);


%% Loop through data
for j = 1:numel(sensors)
    sens = sensors{j};
    R_dir_s = [R_dir, sens, '/'];
    
    for i = Subj
        s_idx = find(i == Subj);
        subj = sprintf('%03d', i);
        
        for k = 1:numel(bands)
            f = bands{k};
            curr = [];
            % get subgraph data
            subset = readNPY([data_dir, subj, '/', sens, '/wpli_', f, '_subset.npy']);
            coeff = readNPY([data_dir, subj, '/',sens, '/wpli_', f, '_coeff.npy']);
            b_exp = subset(:,end);
            [~,bSG] = max(b_exp);
            [~,nbSG] = min(b_exp);
            nSG = size(subset,1);
             
            for n = 1:nSG
                % get expressin into matrix
                curr_SG = subset(n,1:end-1);
                
                flag = min(curr_SG) ~= 0;
                curr = [curr; flag];
            end
            noise_sg{k,i} = curr;
        end        
    end
end

save([save_dir, 'noise_sg.mat'], 'noise_sg')

%% Loop through data - uniform PR

% initialize
noise_sg = cell(numel(bands),nSubj);

for j = 1:numel(sensors)
    sens = sensors{j};
    R_dir_s = [R_dir, sens, '/'];
    
    for i = Subj
        s_idx = find(i == Subj);
        subj = sprintf('%03d', i);
        
        for k = 1:numel(bands)
            f = bands{k};
            curr = [];
            % get subgraph data
            subset = readNPY([data_dir, subj, '/', sens, '/wpli_pr_', f, '_subset.npy']);
            coeff = readNPY([data_dir, subj, '/',sens, '/wpli_pr_', f, '_coeff.npy']);
            b_exp = subset(:,end);
            [~,bSG] = max(b_exp);
            [~,nbSG] = min(b_exp);
            nSG = size(subset,1);
             
            for n = 1:nSG
                % get expressin into matrix
                curr_SG = subset(n,1:end-1);
                
                flag = min(curr_SG) ~= 0;
                curr = [curr; flag];
            end
            noise_sg{k,i} = curr;
        end        
    end
end

save([save_dir, 'pr_noise_sg.mat'], 'noise_sg')

%% Loop through data ind pr

% initialize
noise_sg = cell(numel(bands),nSubj);

for j = 1:numel(sensors)
    sens = sensors{j};
    R_dir_s = [R_dir, sens, '/'];
    
    for i = Subj
        s_idx = find(i == Subj);
        subj = sprintf('%03d', i);
        
        for k = 1:numel(bands)
            f = bands{k};
            curr = [];
            % get subgraph data
            subset = readNPY([data_dir, subj, '/', sens, '/wpli_ind_', f, '_subset.npy']);
            coeff = readNPY([data_dir, subj, '/',sens, '/wpli_ind_', f, '_coeff.npy']);
            b_exp = subset(:,end);
            [~,bSG] = max(b_exp);
            [~,nbSG] = min(b_exp);
            nSG = size(subset,1);
             
            for n = 1:nSG
                % get expressin into matrix
                curr_SG = subset(n,1:end-1);
                
                flag = min(curr_SG) ~= 0;
                curr = [curr; flag];
            end
            noise_sg{k,i} = curr;
        end        
    end
end

save([save_dir, 'ind_noise_sg.mat'], 'noise_sg')