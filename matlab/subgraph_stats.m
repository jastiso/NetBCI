%% Subgraph Stats
% This formats data for exploratory data analysis in R. We will make one
% data structure with the BC, temporal energy, coeff median, subj, band,
% condition and regional expression across all lobes for every subj:subgraph:band
% combo.

addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
R_dir = '/Users/stiso/Documents/R/NetBCI/data/wpli/';

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
% make directories
if ~exist(R_dir, 'dir')
    mkdir(R_dir)
end

Subj = [1:20];
nSubj = numel(Subj);
freqs = [7,14;15,30;31,45;55,70];
bands = [{'alpha'}, {'beta'}, {'low_gamma'}];
sensors = [{'grad'}];
%regions = [{'Left_frontal'}, {'Left_occipital'}, {'Left_parietal'}, {'Left_temporal'}, ...
%    {'Right_frontal'}, {'Right_occipital'}, {'Right_parietal'}, {'Right_temporal'}, {'Vertex'}];
regions = [{'Left_frontal'}, {'Left_occipital'}, {'Left_motor'}, {'Right_frontal'}, {'Right_occipital'}, {'Right_motor'}];
    
nNode = 102;
nEdges = (nNode^2-nNode)/2;

% load behavior
load([top_dir, 'Behavior/stats'])


%% Loop through data

load([save_dir, 'noise_sg.mat']);
noise_idx = noise_sg;

% initialize
sg_n = [];
%noise_sg = []; changed code to  not include this SG
bc = [];
slope = [];
imp = [];
diff = [];
temp_energy = [];
coeff_med = [];
coeff_skew = [];
coeff_min = [];
coeff_quant = [];
coeff_autocorr = [];
subjects = [];
band = [];
cond = [];
Left_frontal = [];
Left_occipital = [];
Left_motor = [];
Right_frontal = [];
Right_motor = [];
Right_occipital = [];

% Left_frontal_wi = [];
% Left_parietal_wi = [];
% Left_temporal_wi = [];
% Left_occipital_wi = [];
% Vertex_wi = [];
% Right_frontal_wi = [];
% Right_parietal_wi = [];
% Right_temporal_wi = [];
% Right_occipital_wi = [];
% Left_frontal_bw = [];
% Left_parietal_bw = [];
% Left_temporal_bw = [];
% Left_occipital_bw = [];
% Vertex_bw = [];
% Right_frontal_bw = [];
% Right_parietal_bw = [];
% Right_temporal_bw = [];
% Right_occipital_bw = [];
exp_max = [];


for j = 1:numel(sensors)
    sens = sensors{j};
    R_dir_s = [R_dir, sens, '/'];
    
    for i = Subj
        s_idx = find(i == Subj);
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
            idx = noise_sg{k,i};
            coeff = coeff(~idx,:);
            subset = subset(~idx,:);
            
            b_exp = subset(:,end);
            nSG = numel(b_exp);
            [tmp,tmp_idx] = sort(b_exp,'descend');
            bSG = tmp_idx(1);
            bSG2 = tmp_idx(2); % second biggest element
            bSG3 = tmp_idx(3); % third biggest element
            [~, nzSG] = min(nonzeros((b_exp))); % smallest nonzero - this is how we will opporationalize "low"
            
            % if you actualy want to look at 0s
            nbSG = find(b_exp == 0);

            % get indices of all relevant subraphs
            %relevant_sgs = [bSG, bSG2, bSG3, nzSG, nbSG];
            %nSG = numel(relevant_sg);
            
            for n = 1:nSG
                % get expressin into matrix
                curr_SG = subset(n,1:end-1);
                node_exp = get_sg_matrix(nNode, curr_SG);
                
                % record behavioral coefficient, subject, and band
                % some subgraphs are both the smallest non-zero, and the
                % thrid highest. how do I deal with this? Right now, they
                % are classified as smallest nonzero...
                
                bc = [bc; subset(n,end)];
                subjects = [subjects; {subj}];
                slope = [slope; betas(i)];
                imp = [imp; improvement(i)];
                diff = [diff; difference(i)];
                band = [band; {f}];
                sg_n = [sg_n; n];
                if n == bSG
                    cond = [cond; {'high'}];
                elseif n == nzSG
                    cond = [cond; {'low'}];
                elseif n == bSG2
                    cond = [cond; {'high2'}];
                elseif n == bSG3
                    cond = [cond; {'high3'}];
                elseif any(find((n == nbSG)))
                    cond = [cond; {'zero'}];
                else
                    cond = [cond; {'other'}];
                end
                
                % temporal stats
                energy = sum(coeff(n,:).^2);
                temp_energy = [temp_energy; energy];
                
                % coeff stats
                coeff_med = [coeff_med; median(curr_SG)];
                coeff_skew = [coeff_skew; skewness(curr_SG)];
                coeff_min = [coeff_min; min(curr_SG)];
                coeff_quant = [coeff_quant; quantile(curr_SG,.025)];
                tmp = morans_I(node_exp, ones(5));
                %pause(0.1)
                coeff_autocorr = [coeff_autocorr; nanmean(nanmean(tmp))];
                
                % regional stats
                for m = 1:numel(regions)
                    load([top_dir, 'montages/', regions{m}, '_idx.mat'])
                    idx = logical(idx);
                    % simplified this to just get the strength. How strong
                    % are thse nodes/all nodes
                    str = (mean(mean(node_exp(idx,:))))/mean(mean(node_exp));
                    eval([regions{m}, ' = [', regions{m}, '; str];']);
                    
                    % distinguishing between withing and between lobe edges
%                     % get the mean edge within a lobe (total edges/number of edges), divided by
%                     % the total edges
%                     mean_within = (sum(sum(node_exp(idx,idx)))/(sum(idx)^2))/sum(sum(node_exp));
%                     % get the mean edge between a lobes (total edges/number of edges), divided by
%                     % the total edges
%                     mean_between = (sum(sum(node_exp(idx,~idx)))/((sum(~idx)*sum(idx))))/sum(sum(node_exp));
%                     eval([regions{m}, '_wi = [', regions{m}, '_wi; mean_within];']);
%                     eval([regions{m}, '_bw = [', regions{m}, '_bw; mean_between];']);
                end
                % specifically frontal to vertex
%                 load([top_dir, 'montages/Right_frontal_idx.mat'])
%                 rf_idx = idx;
%                 load([top_dir, 'montages/Left_frontal_idx.mat'])
%                 lf_idx = idx;
%                 load([top_dir, 'montages/Vertex_idx.mat'])
%                 v_idx = idx;
%                 r_v = (sum(sum(node_exp(rf_idx,v_idx)))/((sum(rf_idx)*sum(v_idx))))/sum(sum(node_exp));
%                 l_v = (sum(sum(node_exp(lf_idx,v_idx)))/((sum(lf_idx)*sum(v_idx))))/sum(sum(node_exp));
%                 rfv = [rfv; r_v];
%                 lfv = [lfv; l_v];
                
            end
        end
    end
end

%save([R_dir, 'subgraph_stats.mat'], 'slope', 'imp', 'diff', 'noise_sg', 'sg_n', 'subjects', 'band', 'cond', 'bc', 'temp_energy', ...
%    'coeff_med', 'coeff_skew', 'coeff_min', 'coeff_autocorr', 'coeff_quant', 'Left_frontal_wi', 'Left_occipital_wi',...
%    'Left_parietal_wi', 'Left_temporal_wi', 'Right_frontal_wi', 'Right_occipital_wi', 'Right_parietal_wi', 'Right_temporal_wi', 'Vertex_wi', ...
%    'Left_frontal_bw', 'Left_occipital_bw','Left_parietal_bw', 'Left_temporal_bw', 'Right_frontal_bw', 'Right_occipital_bw', 'Right_parietal_bw',...
%    'Right_temporal_bw', 'Vertex_bw', 'rfv', 'lfv', 'exp_max')

save([R_dir, 'subgraph_stats.mat'], 'slope', 'imp', 'diff', 'noise_sg', 'sg_n', 'subjects', 'band', 'cond', 'bc', 'temp_energy', ...
    'coeff_med', 'coeff_skew', 'coeff_min', 'coeff_autocorr', 'coeff_quant', 'Left_frontal', 'Left_occipital','Left_motor', ...
    'Right_frontal', 'Right_occipital', 'Right_motor', 'exp_max')


%% Loop through data - uniform pr

load([save_dir, 'pr_noise_sg.mat']);
noise_idx = noise_sg;

% initialize
sg_n = [];
%noise_sg = []; changed code to  not include this SG
bc = [];
slope = [];
imp = [];
diff = [];
temp_energy = [];
coeff_med = [];
coeff_skew = [];
coeff_min = [];
coeff_quant = [];
coeff_autocorr = [];
subjects = [];
band = [];
cond = [];
Left_frontal = [];
Left_occipital = [];
Left_motor = [];
Right_frontal = [];
Right_motor = [];
Right_occipital = [];

% Left_frontal_wi = [];
% Left_parietal_wi = [];
% Left_temporal_wi = [];
% Left_occipital_wi = [];
% Vertex_wi = [];
% Right_frontal_wi = [];
% Right_parietal_wi = [];
% Right_temporal_wi = [];
% Right_occipital_wi = [];
% Left_frontal_bw = [];
% Left_parietal_bw = [];
% Left_temporal_bw = [];
% Left_occipital_bw = [];
% Vertex_bw = [];
% Right_frontal_bw = [];
% Right_parietal_bw = [];
% Right_temporal_bw = [];
% Right_occipital_bw = [];
exp_max = [];


for j = 1:numel(sensors)
    sens = sensors{j};
    R_dir_s = [R_dir, sens, '/'];
    
    for i = Subj
        s_idx = find(i == Subj);
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
            idx = noise_sg{k,i};
            coeff = coeff(~idx,:);
            subset = subset(~idx,:);
            nSG = size(coeff,1);

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
            
            for n = 1:nSG
                % get expressin into matrix
                curr_SG = subset(n,1:end-1);
                node_exp = get_sg_matrix(nNode, curr_SG);
                
                % record behavioral coefficient, subject, and band
                bc = [bc; subset(n,end)];
                subjects = [subjects; {subj}];
                slope = [slope; betas(i)];
                imp = [imp; improvement(i)];
                diff = [diff; difference(i)];
                band = [band; {f}];
                sg_n = [sg_n; n];
                if n == bSG
                    cond = [cond; {'high'}];
                elseif n == nzSG
                    cond = [cond; {'low'}];
                elseif n == bSG2
                    cond = [cond; {'high2'}];
                elseif n == bSG3
                    cond = [cond; {'high3'}];
                elseif any(find((n == nbSG)))
                    cond = [cond; {'zero'}];
                else
                    cond = [cond; {'other'}];
                end
                
                % temporal stats
                energy = sum(coeff(n,:).^2);
                temp_energy = [temp_energy; energy];
                
                % coeff stats
                coeff_med = [coeff_med; median(curr_SG)];
                coeff_skew = [coeff_skew; skewness(curr_SG)];
                coeff_min = [coeff_min; min(curr_SG)];
                coeff_quant = [coeff_quant; quantile(curr_SG,.025)];
                tmp = morans_I(node_exp, ones(5));
                %pause(0.1)
                coeff_autocorr = [coeff_autocorr; nanmean(nanmean(tmp))];
                
                % regional stats
                for m = 1:numel(regions)
                    load([top_dir, 'montages/', regions{m}, '_idx.mat'])
                    idx = logical(idx);
                    % simplified this to just get the strength. How strong
                    % are thse nodes/all nodes
                    str = (mean(mean(node_exp(idx,:))))/mean(mean(node_exp));
                    eval([regions{m}, ' = [', regions{m}, '; str];']);
                    
                    % distinguishing between withing and between lobe edges
%                     % get the mean edge within a lobe (total edges/number of edges), divided by
%                     % the total edges
%                     mean_within = (sum(sum(node_exp(idx,idx)))/(sum(idx)^2))/sum(sum(node_exp));
%                     % get the mean edge between a lobes (total edges/number of edges), divided by
%                     % the total edges
%                     mean_between = (sum(sum(node_exp(idx,~idx)))/((sum(~idx)*sum(idx))))/sum(sum(node_exp));
%                     eval([regions{m}, '_wi = [', regions{m}, '_wi; mean_within];']);
%                     eval([regions{m}, '_bw = [', regions{m}, '_bw; mean_between];']);
                end
                % specifically frontal to vertex
%                 load([top_dir, 'montages/Right_frontal_idx.mat'])
%                 rf_idx = idx;
%                 load([top_dir, 'montages/Left_frontal_idx.mat'])
%                 lf_idx = idx;
%                 load([top_dir, 'montages/Vertex_idx.mat'])
%                 v_idx = idx;
%                 r_v = (sum(sum(node_exp(rf_idx,v_idx)))/((sum(rf_idx)*sum(v_idx))))/sum(sum(node_exp));
%                 l_v = (sum(sum(node_exp(lf_idx,v_idx)))/((sum(lf_idx)*sum(v_idx))))/sum(sum(node_exp));
%                 rfv = [rfv; r_v];
%                 lfv = [lfv; l_v];
                
            end
        end
    end
end

%save([R_dir, 'subgraph_stats.mat'], 'slope', 'imp', 'diff', 'noise_sg', 'sg_n', 'subjects', 'band', 'cond', 'bc', 'temp_energy', ...
%    'coeff_med', 'coeff_skew', 'coeff_min', 'coeff_autocorr', 'coeff_quant', 'Left_frontal_wi', 'Left_occipital_wi',...
%    'Left_parietal_wi', 'Left_temporal_wi', 'Right_frontal_wi', 'Right_occipital_wi', 'Right_parietal_wi', 'Right_temporal_wi', 'Vertex_wi', ...
%    'Left_frontal_bw', 'Left_occipital_bw','Left_parietal_bw', 'Left_temporal_bw', 'Right_frontal_bw', 'Right_occipital_bw', 'Right_parietal_bw',...
%    'Right_temporal_bw', 'Vertex_bw', 'rfv', 'lfv', 'exp_max')

save([R_dir, 'pr_subgraph_stats.mat'], 'slope', 'imp', 'diff', 'noise_sg', 'sg_n', 'subjects', 'band', 'cond', 'bc', 'temp_energy', ...
    'coeff_med', 'coeff_skew', 'coeff_min', 'coeff_autocorr', 'coeff_quant', 'Left_frontal', 'Left_occipital','Left_motor', ...
    'Right_frontal', 'Right_occipital', 'Right_motor', 'exp_max')