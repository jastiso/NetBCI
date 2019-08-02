%% Get correlation between expression and behavior
% mostly a sonity check

addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath(genpath('/Users/stiso/Documents/MATLAB/mult_comp_perm_corr/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
R_dir_o = '/Users/stiso/Documents/R/NetBCI/data/wpli/';

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/param/';
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

Subj = [1:20];
nSubj = numel(Subj);
freqs = [7,14;15,30;31,45;55,70];
bands = [{'alpha'},{'beta'},{'low_gamma'}];

sensors = [{'grad'}];
% load behavior
load([top_dir, 'Behavior/behavior_all'])
load([top_dir, 'Behavior/stats'])
load([save_dir, 'noise_sg.mat']);

%initialize
max_exp = zeros(nSubj, numel(bands), numel(sensors));
num_zero = zeros(nSubj, numel(bands), numel(sensors));
max_exp2 = zeros(nSubj, numel(bands), numel(sensors));
max_exp3 = zeros(nSubj, numel(bands), numel(sensors));
max_exp4 = zeros(nSubj, numel(bands), numel(sensors));
min_exp = zeros(nSubj, numel(bands), numel(sensors));
min_exp2 = zeros(nSubj, numel(bands), numel(sensors));
sum_exp = zeros(nSubj, numel(bands), numel(sensors));
sd_exp = zeros(nSubj, numel(bands), numel(sensors));

%% Loop through data

for j = 1:numel(sensors)
    sens = sensors{j};
    R_dir = [R_dir_o, sens, '/'];
    
    for i = Subj
        idx = find(i == Subj);
        subj = sprintf('%03d', i);
        
        for k = 1:numel(bands)
            f = bands{k};
            
            % get subgraph data
            if strcmp(data_dir(end-5:end-1),'param')
                subset = readNPY([data_dir, subj, '/', sens, 'wpli_', f, '_subset.npy']);
                coeff = readNPY([data_dir, subj, '/',sens, 'wpli_', f, '_coeff.npy']);
            else
                subset = readNPY([data_dir, subj, '/', sens, '/wpli_', f, '_subset.npy']);
                coeff = readNPY([data_dir, subj, '/',sens, '/wpli_', f, '_coeff.npy']);
            end
            
            % remove noise SG
            nidx = noise_sg{k,i};
            coeff = coeff(~nidx,:);
            subset = subset(~nidx,:);
            
            %get subgraph with highest expression
            [max_e,bsg] = max(subset(:,end));
            b_exp = subset(:,end);
            [max_e2,bsg2] = max(b_exp(b_exp<max(b_exp))); % second biggest element
            [tmp,tmp_idx] = sort(b_exp,'descend');
            max_e3 = tmp(3); % third biggest element
            bsg3 = tmp_idx(3); % third biggest element
            max_e4 = tmp(4); % 4th biggest element
            bsg4 = tmp_idx(4); % 4th biggest element
            [min_e,nbsg] = min(nonzeros(subset(:,end))); % here, low is smallest nonzero
            [tmp,tmp_idx] = sort(nonzeros(subset(:,end)),'ascend');
            try
                min_e2 = tmp(2); % second smallest element
                nbsg2 = tmp_idx(2); % second smallest  element
            catch
                min_e2 = tmp(1); % if only one non_zero
                nbsg2 = tmp_idx(1);
            end
            
            % add to structure
            num_zero(idx,k,j) = sum(b_exp == 0);
            max_exp(idx,k,j) = max_e;
            max_exp2(idx,k,j) = max_e2;
            max_exp3(idx,k,j) = max_e3;
            max_exp4(idx,k,j) = max_e4;
            min_exp(idx,k,j) = min_e;
            min_exp2(idx,k,j) = min_e2;
            sum_exp(idx,k,j) = mean(nonzeros(subset(:,end)));
            sd_exp(idx,k,j) = std(nonzeros(subset(:,end)));
            
            % visual check that there is some distribution of performance
            % loadings
            figure(1); clf
            hist(subset(:,end))
            pause(0.001)
            %saveas(gca, [save_dir, subj, '_', f, '_exp_hist.png'], 'png')
        end
    end
    
    % save to R_dir
    save([R_dir, 'exp_beh_cor.mat'], 'Subj', 'betas', 'num_zero', 'max_exp', 'max_exp2', 'min_exp', 'min_exp2', 'max_exp3', 'max_exp4', 'sum_exp', 'sd_exp');
    
end


%% Null models - uniform PR

%load noise sg
load([save_dir, 'pr_noise_sg.mat']);

%initialize
max_exp_null = zeros(nSubj, numel(bands), numel(sensors));
max_exp_null2 = zeros(nSubj, numel(bands), numel(sensors));
max_exp_null3 = zeros(nSubj, numel(bands), numel(sensors));
max_exp_null4 = zeros(nSubj, numel(bands), numel(sensors));
min_exp_null = zeros(nSubj, numel(bands), numel(sensors));
min_exp_null2 = zeros(nSubj, numel(bands), numel(sensors));
sum_exp_null = zeros(nSubj, numel(bands), numel(sensors));
sd_exp_null = zeros(nSubj, numel(bands), numel(sensors));


% Loop through data
for j = 1:numel(sensors)
    sens = sensors{j};
    R_dir = [R_dir_o, sens, '/'];
    
    for i = Subj
        idx = find(i == Subj);
        subj = sprintf('%03d', i);
        
        for k = 1:numel(bands)
            f = bands{k};
            
            % get subgraph data
            if strcmp(data_dir(end-5:end-1),'param')
                subset = readNPY([data_dir, subj, '/', sens, 'wpli_pr_', f, '_subset.npy']);
                coeff = readNPY([data_dir, subj, '/',sens, 'wpli_pr_', f, '_coeff.npy']);
            else
                subset = readNPY([data_dir, subj, '/', sens, '/wpli_pr_', f, '_subset.npy']);
                coeff = readNPY([data_dir, subj, '/',sens, '/wpli_pr_', f, '_coeff.npy']);
            end
            
            % remove noise SG
            nidx = noise_sg{k,i};
            coeff = coeff(~nidx,:);
            subset = subset(~nidx,:);
            
            %get subgraph with highest expression
            [max_e,bsg] = max(subset(:,end));
            b_exp = subset(:,end);
            [max_e2,bsg2] = max(b_exp(b_exp<max(b_exp))); % second biggest element
            [tmp,tmp_idx] = sort(b_exp,'descend');
            max_e3 = tmp(3); % third biggest element
            bsg3 = tmp_idx(3); % third biggest element
            max_e4 = tmp(4); % 4th biggest element
            bsg4 = tmp_idx(4); % 4th biggest element
            [min_e,nbsg] = min(nonzeros(subset(:,end))); % here, low is smallest nonzero
            [tmp,tmp_idx] = sort(nonzeros(subset(:,end)),'ascend');
            try
                min_e2 = tmp(2); % second smallest element
                nbsg2 = tmp_idx(2); % second smallest  element
            catch
                min_e2 = tmp(1); % if only one non_zero
                nbsg2 = tmp_idx(1);
            end
            
            % get corr
            max_exp_null(idx,k,j) = max_e;
            max_exp_null2(idx,k,j) = max_e2;
            max_exp_null3(idx,k,j) = max_e3;
            max_exp_null4(idx,k,j) = max_e4;
            min_exp_null(idx,k,j) = min_e;
            min_exp_null2(idx,k,j) = min_e2;
            sum_exp_null(idx,k,j) = mean(nonzeros(subset(:,end)));
            sd_exp_null(idx,k,j) = std(nonzeros(subset(:,end)));
        end
    end

    % save to R_dir
    save([R_dir, 'exp_beh_cor_pr.mat'], 'Subj', 'betas',  'max_exp_null', 'max_exp_null2', 'max_exp_null3', 'max_exp_null4', 'min_exp_null', 'min_exp_null2', 'sum_exp_null', 'sd_exp_null');
    
end

