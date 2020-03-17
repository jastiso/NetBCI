%% Compare resting state FC to each subgraph

clear
clc

addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
addpath('/Users/stiso/Documents/MATLAB/control_example/')
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
freqs = [7,14;15,30;31,45];
bands = [{'alpha'}, {'beta'}, {'low_gamma'}];
sensors = [{'grad'}];


nNode = 102;
nEdges = (nNode^2-nNode)/2;
load([save_dir, 'noise_sg.mat']);

% initialize
rs_sim = cell2table(cell(0,4), 'VariableNames', {'sim', 'subj', 'band', 'sg'});

%% Loop through data

sens = 'grad';

cnt = 1; % for order
for i = Subj
    s_idx = find(i == Subj);
    subj = sprintf('%03d', i);
    
    for k = 1:numel(bands)
        f = bands{k};
        fprintf('\nSubj %03d... band %s...', i, f)

        % get subgraph data
        if strcmp(data_dir(end-5:end-1),'param')
            subset = readNPY([data_dir, subj, '/', sens, 'wpli_', f, '_subset.npy']);
            coeff = readNPY([data_dir, subj, '/',sens, 'wpli_', f, '_coeff.npy']);
        else
            subset = readNPY([data_dir, subj, '/', sens, '/wpli_', f, '_subset.npy']);
            coeff = readNPY([data_dir, subj, '/',sens, '/wpli_', f, '_coeff.npy']);
        end
        
        % remove noise SG
        idx = noise_sg{k,i};
        coeff = coeff(~idx,:);
        subset = subset(~idx,:);
        
        b_exp = subset(:,end);
        [tmp,tmp_idx] = sort(b_exp,'descend');
        bSG = tmp_idx(1);
        bSG2 = tmp_idx(2); % second biggest element
        bSG3 = tmp_idx(3); % third biggest element
        [~, nzSG] = min(nonzeros((b_exp))); % smallest nonzero - this is how we will opporationalize "low"
        
        % if you actualy want to look at 0s
        if sum(b_exp == 0) <= 1
            [~,nbSG] = min(b_exp);
        else
            nbSG = find(b_exp == 0);
        end
        nZero = numel(nbSG);
        
        % optimal control
        % scale to be stable
        % 1st
        high_mat = subset(bSG,(1:end-1));
        %2nd
        high2_mat = subset(bSG2,1:(end-1));
        %3rd
        high3_mat = subset(bSG3,1:(end-1));
        %lowest
        low_mat = subset(nzSG,1:(end-1));
        
        % get rs
        load([top_dir, 'Session4/RS2_END/', subj, '/FCmatrices/', f, '_rs_wpli.mat']);
        wpli = wpli(logical(triu(ones(nNode),1)));
        
        % cosine similarity with rs
        high_sim = 1 - pdist([high_mat(high_mat ~= 0); wpli(high_mat ~= 0)'],'cosine');
        high2_sim = 1 - pdist([high2_mat(high2_mat ~= 0); wpli(high2_mat ~= 0)'],'cosine');
        high3_sim = 1 - pdist([high3_mat(high3_mat ~= 0); wpli(high3_mat ~= 0)'],'cosine');
        low_sim = 1 - pdist([low_mat(low_mat ~= 0); wpli(low_mat ~= 0)'],'cosine');
        
        % add to table
        rs_sim = [rs_sim; table([high_sim; high2_sim; high3_sim; low_sim], repmat(subj,4,1),...
            repmat({f}, 4,1), {'high'; 'high2'; 'high3'; 'low'}, 'VariableNames', {'sim', 'subj', 'band', 'sg'})];
    end
end

writetable(rs_sim, [R_dir, 'grad/rs_sim.csv'])

