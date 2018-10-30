%% Compare parameters
% Check if selected parameters are stable across runs of NMF

addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
R_dir = '/Users/stiso/Documents/R/NetBCI/data/wpli/';

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir1 = '/Users/stiso/Documents/Python/NetBCI/NMF/param/';
data_dir2 = '/Users/stiso/Documents/Python/NetBCI/NMF/';
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
bands = [{'alpha'}, {'beta'}, {'low_gamma'}];
sensors = [{'grad'}];

%% Loop through subjects

alpha = zeros(nSubj, 2, numel(bands));
beta = zeros(nSubj, 2, numel(bands));
m = zeros(nSubj, 2, numel(bands));

sens = sensors{1};
for i = Subj
    s_idx = find(i == Subj);
    subj = sprintf('%03d', i);
    
    for k = 1:numel(bands)
        f = bands{k};
        
        % get subgraph data
        param1 = readNPY([data_dir1, subj, '/', sens, 'wpli_', f, '_params.npy']);
        param2 = readNPY([data_dir2, subj, '/',sens, '/wpli_', f, '_params.npy']);
        
        
    end
end